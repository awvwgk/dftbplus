!> Add-on for tblite to enable lambda scaling calculations
module dftbp_extlibs_tblite_lambda
  use dftbp_common_accuracy, only : dp
  use mctc_io, only : structure_type
  use multicharge, only : get_coordination_number, get_lattice_points
  use multicharge_data, only : get_covalent_rad
  implicit none
  private

  public :: get_fragment, unfold_fragment, guess_bonds

  integer, parameter :: max_elem = 86
  real(dp), parameter :: en(max_elem) = [&
    & 2.30085633_dp, 2.78445145_dp, 1.52956084_dp, 1.51714704_dp, 2.20568300_dp, &
    & 2.49640820_dp, 2.81007174_dp, 4.51078438_dp, 4.67476223_dp, 3.29383610_dp, &
    & 2.84505365_dp, 2.20047950_dp, 2.31739628_dp, 2.03636974_dp, 1.97558064_dp, &
    & 2.13446570_dp, 2.91638164_dp, 1.54098156_dp, 2.91656301_dp, 2.26312147_dp, &
    & 2.25621439_dp, 1.32628677_dp, 2.27050569_dp, 1.86790977_dp, 2.44759456_dp, &
    & 2.49480042_dp, 2.91545568_dp, 3.25897750_dp, 2.68723778_dp, 1.86132251_dp, &
    & 2.01200832_dp, 1.97030722_dp, 1.95495427_dp, 2.68920990_dp, 2.84503857_dp, &
    & 2.61591858_dp, 2.64188286_dp, 2.28442252_dp, 1.33011187_dp, 1.19809388_dp, &
    & 1.89181390_dp, 2.40186898_dp, 1.89282464_dp, 3.09963488_dp, 2.50677823_dp, &
    & 2.61196704_dp, 2.09943450_dp, 2.66930105_dp, 1.78349472_dp, 2.09634533_dp, &
    & 2.00028974_dp, 1.99869908_dp, 2.59072029_dp, 2.54497829_dp, 2.52387890_dp, &
    & 2.30204667_dp, 1.60119300_dp, 2.00000000_dp, 2.00000000_dp, 2.00000000_dp, &
    & 2.00000000_dp, 2.00000000_dp, 2.00000000_dp, 2.00000000_dp, 2.00000000_dp, &
    & 2.00000000_dp, 2.00000000_dp, 2.00000000_dp, 2.00000000_dp, 2.00000000_dp, &
    & 2.00000000_dp, 2.30089349_dp, 1.75039077_dp, 1.51785130_dp, 2.62972945_dp, &
    & 2.75372921_dp, 2.62540906_dp, 2.55860939_dp, 3.32492356_dp, 2.65140898_dp, &
    & 1.52014458_dp, 2.54984804_dp, 1.72021963_dp, 2.69303422_dp, 1.81031095_dp, &
    & 2.34224386_dp]
  real(dp), parameter :: r0(max_elem) = [&
    & 0.55682207_dp, 0.80966997_dp, 2.49092101_dp, 1.91705642_dp, 1.35974851_dp, &
    & 0.98310699_dp, 0.98423007_dp, 0.76716063_dp, 1.06139799_dp, 1.17736822_dp, &
    & 2.85570926_dp, 2.56149012_dp, 2.31673425_dp, 2.03181740_dp, 1.82568535_dp, &
    & 1.73685958_dp, 1.97498207_dp, 2.00136196_dp, 3.58772537_dp, 2.68096221_dp, &
    & 2.23355957_dp, 2.33135502_dp, 2.15870365_dp, 2.10522128_dp, 2.16376162_dp, &
    & 2.10804037_dp, 1.96460045_dp, 2.00476257_dp, 2.22628712_dp, 2.43846700_dp, &
    & 2.39408483_dp, 2.24245792_dp, 2.05751204_dp, 2.15427677_dp, 2.27191920_dp, &
    & 2.19722638_dp, 3.80910350_dp, 3.26020971_dp, 2.99716916_dp, 2.71707818_dp, &
    & 2.34950167_dp, 2.11644818_dp, 2.47180659_dp, 2.32198800_dp, 2.32809515_dp, &
    & 2.15244869_dp, 2.55958313_dp, 2.59141300_dp, 2.62030465_dp, 2.39935278_dp, &
    & 2.56912355_dp, 2.54374096_dp, 2.56914830_dp, 2.53680807_dp, 4.24537037_dp, &
    & 3.66542289_dp, 3.19903011_dp, 2.80000000_dp, 2.80000000_dp, 2.80000000_dp, &
    & 2.80000000_dp, 2.80000000_dp, 2.80000000_dp, 2.80000000_dp, 2.80000000_dp, &
    & 2.80000000_dp, 2.80000000_dp, 2.80000000_dp, 2.80000000_dp, 2.80000000_dp, &
    & 2.80000000_dp, 2.34880037_dp, 2.37597108_dp, 2.49067697_dp, 2.14100577_dp, &
    & 2.33473532_dp, 2.19498900_dp, 2.12678348_dp, 2.34895048_dp, 2.33422774_dp, &
    & 2.86560827_dp, 2.62488837_dp, 2.88376127_dp, 2.75174124_dp, 2.83054552_dp, &
    & 2.63264944_dp]
  real(dp), parameter :: cnfak(max_elem) = [&
    & 0.17957827_dp, 0.25584045_dp,-0.02485871_dp, 0.00374217_dp, 0.05646607_dp, &
    & 0.10514203_dp, 0.09753494_dp, 0.30470380_dp, 0.23261783_dp, 0.36752208_dp, &
    & 0.00131819_dp,-0.00368122_dp,-0.01364510_dp, 0.04265789_dp, 0.07583916_dp, &
    & 0.08973207_dp,-0.00589677_dp, 0.13689929_dp,-0.01861307_dp, 0.11061699_dp, &
    & 0.10201137_dp, 0.05426229_dp, 0.06014681_dp, 0.05667719_dp, 0.02992924_dp, &
    & 0.03764312_dp, 0.06140790_dp, 0.08563465_dp, 0.03707679_dp, 0.03053526_dp, &
    &-0.00843454_dp, 0.01887497_dp, 0.06876354_dp, 0.01370795_dp,-0.01129196_dp, &
    & 0.07226529_dp, 0.01005367_dp, 0.01541506_dp, 0.05301365_dp, 0.07066571_dp, &
    & 0.07637611_dp, 0.07873977_dp, 0.02997732_dp, 0.04745400_dp, 0.04582912_dp, &
    & 0.10557321_dp, 0.02167468_dp, 0.05463616_dp, 0.05370913_dp, 0.05985441_dp, &
    & 0.02793994_dp, 0.02922983_dp, 0.02220438_dp, 0.03340460_dp,-0.04110969_dp, &
    &-0.01987240_dp, 0.07260201_dp, 0.07700000_dp, 0.07700000_dp, 0.07700000_dp, &
    & 0.07700000_dp, 0.07700000_dp, 0.07700000_dp, 0.07700000_dp, 0.07700000_dp, &
    & 0.07700000_dp, 0.07700000_dp, 0.07700000_dp, 0.07700000_dp, 0.07700000_dp, &
    & 0.07700000_dp, 0.08379100_dp, 0.07314553_dp, 0.05318438_dp, 0.06799334_dp, &
    & 0.04671159_dp, 0.06758819_dp, 0.09488437_dp, 0.07556405_dp, 0.13384502_dp, &
    & 0.03203572_dp, 0.04235009_dp, 0.03153769_dp,-0.00152488_dp, 0.02714675_dp, &
    & 0.04800662_dp]

contains

  subroutine guess_bonds(mol, bond)
    type(structure_type), intent(in) :: mol
    real(dp), intent(out) :: bond(:, :)

    integer, allocatable :: num(:)
    real(dp), allocatable :: rcov(:), cn(:), r0(:, :), lattr(:, :)
    real(dp), parameter :: cutoff = 20.0_dp
    integer :: iat, jat, itr
    real(dp) :: r1, m1, vec(3)

    allocate(r0(mol%nat, mol%nat), cn(mol%nat))

    rcov = get_covalent_rad(mol%num)
    call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
    call get_coordination_number(mol, lattr, cutoff, rcov, cn=cn)

    num = mol%num(mol%id)
    call gfnffrab(mol%nat, num, cn, r0)

    do iat = 1, mol%nat
      do jat = 1, iat - 1
        m1 = huge(1.0_dp)
        do itr = 1, size(lattr, 2)
          vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - lattr(:, itr)
          m1 = min(norm2(vec), m1)
        end do
        bond(jat, iat) = erf_count(7.5_dp, m1, r0(jat, iat))
        bond(iat, jat) = bond(jat, iat)
      end do
    end do

  end subroutine guess_bonds


  !> Error function counting function for coordination number contributions.
  pure function erf_count(k, r, r0) result(count)

    !> Steepness of the counting function.
    real(dp), intent(in) :: k

    !> Current distance.
    real(dp), intent(in) :: r

    !> Cutoff radius.
    real(dp), intent(in) :: r0

    real(dp) :: count

    count = 0.5_dp * (1.0_dp + erf(-k*(r-r0)/r0))

  end function erf_count


  subroutine gfnffrab(n,at,cn,rab)
    !> number of atoms
    integer, intent(in) ::  n
    !> Ordinal number
    integer, intent(in) :: at(:)
    !> Coordination number
    real(dp), intent(in) :: cn(:)
    !> output bond lengths estimates
    real(dp), intent(out) :: rab(:, :)

    integer ::  m,i,j,ii,jj,ati,atj,ir,jr

    real(dp) :: ra,rb,k1,k2,den,ff

    real(dp), parameter :: p(6, 2) = 0.01_dp * reshape([&
      & 29.84522887_dp,-1.70549806_dp, 6.54013762_dp, 6.39169003_dp, 6.00_dp, 5.6_dp, &
      & -8.87843763_dp, 2.10878369_dp, 0.08009374_dp,-0.85808076_dp,-1.15_dp,-1.3_dp],&
      & shape(p))

    do i=1,n
      do j=1,i
        ati=at(i)
        atj=at(j)
        ir=itabrow6(ati)
        jr=itabrow6(atj)
        ra=r0(ati)+cnfak(ati)*cn(i)
        rb=r0(atj)+cnfak(atj)*cn(j)
        den=abs(en(ati)-en(atj))
        k1=0.5_dp*(p(ir,1)+p(jr,1))
        k2=0.5_dp*(p(ir,2)+p(jr,2))
        ff=1.0_dp-k1*den-k2*den**2
        rab(j,i)=(ra+rb)*ff
        rab(i,j)=rab(j,i)
      enddo
    enddo

  contains
    pure function itabrow6(i)
      integer, intent(in) :: i
      integer :: itabrow6

      itabrow6=0
      if (i.gt. 0 .and. i.le. 2) then
        itabrow6=1
      else if (i.gt. 2 .and. i.le.10) then
        itabrow6=2
      else if (i.gt.10 .and. i.le.18) then
        itabrow6=3
      else if (i.gt.18 .and. i.le.36) then
        itabrow6=4
      else if (i.gt.36 .and. i.le.54) then
        itabrow6=5
      else if (i.gt.54) then
        itabrow6=6
      end if
    end function itabrow6
  end subroutine gfnffrab


  subroutine unfold_fragment(abc, wbo, thr)
    !> Bond orders
    real(dp),intent(in) :: wbo(:, :)
    !> Threshold to count bond orders as actual bond
    real(dp),intent(in) :: thr
    !> Fractional coordinates
    real(dp), intent(inout) :: abc(:, :)

    real(dp),allocatable :: bond(:, :)
    logical, allocatable :: taken(:)
    integer :: i, nfrag, nat
    real(dp) :: xsum

    nat = size(abc, 2)
    allocate(taken(nat), source=.false.)
    allocate(bond(nat, nat), source=0.0_dp)
    where(wbo > thr)
      bond = min(wbo, 1.0_dp)
    elsewhere
      bond = 0.0_dp
    endwhere
    nfrag = 0
    do i = 1, nat
      if(taken(i)) cycle
      taken(i)=.true.
      call unfold_neighbours(i, sum(ceiling(bond(:,:)), 1), taken, bond, abc)
    end do
  end subroutine unfold_fragment


  !> Worker routine to transverse the graph recursively
  recursive subroutine unfold_neighbours(i, cn, taken, bond, abc)
    !> Current node in the graph
    integer, intent(in) :: i
    !> Number of edges for every node
    integer, intent(in) :: cn(:)
    !> Status of nodes
    logical, intent(inout) :: taken(:)
    !> Weights of the graph edges, destroyed while transversing the graph
    real(dp), intent(inout) :: bond(:,:)
    !> Fractional coordinates
    real(dp), intent(inout) :: abc(:, :)

    integer :: j, k

    do k = 1, cn(i)
      j = maxloc(bond(:, i), 1)
      bond(j, i) = 0
      if (i == j .or. taken(j)) cycle
      taken(j) = .true.
      abc(:, j) = abc(:, j) - nint(abc(:, j) - abc(:, i))
      call unfold_neighbours(j, cn, taken, bond, abc)
    end do
  end subroutine unfold_neighbours


  !> Obtain fragment information from Wiberg/Mayer type bond orders
  subroutine get_fragment(fragment, wbo, thr)
    !> Fragment information
    integer, intent(out) :: fragment(:)
    !> Bond orders
    real(dp),intent(in) :: wbo(:, :)
    !> Threshold to count bond orders as actual bond
    real(dp),intent(in) :: thr

    real(dp),allocatable :: bond(:, :)
    logical, allocatable :: taken(:)
    integer :: i, nfrag, nat
    real(dp) :: xsum

    nat = size(fragment)
    fragment(:) = 0
    allocate(taken(nat), source=.false.)
    allocate(bond(nat, nat), source=0.0_dp)
    where(wbo > thr)
      bond = min(wbo, 1.0_dp)
    elsewhere
      bond = 0.0_dp
    endwhere

    nfrag = 0
    do i = 1, nat
      if(taken(i)) cycle
      nfrag=nfrag+1
      fragment(i)=nfrag
      taken(i)=.true.
      call find_neighbours(i, sum(ceiling(bond(:,:)), 1), taken, bond, fragment, nfrag)
    end do

  end subroutine get_fragment

  !> Worker routine to transverse the graph recursively
  recursive subroutine find_neighbours(i, cn, taken, bond, fragment, this_fragment)
    !> Current node in the graph
    integer, intent(in) :: i
    !> Number of edges for every node
    integer, intent(in) :: cn(:)
    !> Status of nodes
    logical, intent(inout) :: taken(:)
    !> Weights of the graph edges, destroyed while transversing the graph
    real(dp), intent(inout) :: bond(:,:)
    !> Fragments found
    integer, intent(inout) :: fragment(:)
    !> Index of the current fragmant
    integer, intent(inout) :: this_fragment

    integer :: j, k, dir
    real(dp) :: vec(3)

    do k = 1, cn(i)
      j = maxloc(bond(:, i), 1)
      bond(j, i) = 0
      if (i == j .or. taken(j)) cycle
      fragment(j) = this_fragment
      taken(j) = .true.
      call find_neighbours(j, cn, taken, bond, fragment, this_fragment)
    end do
  end subroutine find_neighbours

end module dftbp_extlibs_tblite_lambda
