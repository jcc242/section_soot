module soot_params
  use Types, only: wp=>dp
  use Bins
  implicit none

  !! These are all read from the input file
  integer :: ND
  !! Number of discrete bins
  integer :: bins_per_2
  !! 2^(1//bins_per_2) is the geometric scalar factor between neighbor sectional bins.
  !! If you do not want any sectional bins, set this to zero
  integer :: n_step
  !! Number of time steps
  real(wp) total_time
  !! Total simulation time
  real(wp) :: w_loss
  !! Constant part of wall loss rate, in s^-1
  real(wp) :: scg_loss
  !! Loss constant of scavenging by pre-existing particles
  real(wp) :: dilution
  !! Dilution constant, in s^-1
  real(wp) :: sat_pressure
  !! Saturation vapor pressure
  logical :: const_monomer
  !! If true, monomer concetration will remain constant
  logical :: use_coag
  !! If true, enable coagulation between particles
  logical :: calc_coefficients
  !! If true, recalculate coefficients between bins, otherwise laod previous coefficients
  logical :: calc_sink_params
  !! If true, recalculate coefficients for sink processes, otherwise load previous coefficients
  logical :: grow_bey_bound
  !! If true, particles are allowed to grow out of the upper boundary
  integer :: kernel
  !! Kernel choice, 1=fuchs, 2=free, 3=hogan, 4=sceats, 5=uniform
  real(wp) :: RR
  !! Monomer production rate, in #/(cm^3*s)
  real(wp), dimension(:), allocatable :: init_conc
  !! Initial number distribution of particles
  real(wp) :: temperature
  !! Temperature in Kelvin
  real(wp) :: pressure
  !! Pressure in Pascal
  real(wp) :: rho
  !! Bulk density in kg/m^3
  real(wp) :: m_mono
  !! Monomer mass in kg
  real(wp) :: surface_tension
  !! Surface tension of nucleation species, in N/m


  !! These are derived quantities
  integer :: DIS_SECT
  !! Selects the sectional model type.
  !! 1 is a discrete sectional model, 2 is a pure sectional model,
  !! and 3 is a pure discrete model.
  integer :: NS
  !! Number of sectional bins.
  !! If you do not want any sectional bins, set bins_per_2 to zero
  integer NT
  !! Number of total bins
  real(wp) :: v_mono
  !! Volume of monomer
  real(wp) :: d_mono
  !! Diameter of monomer
  real(wp) :: A
  !! Dimensionless surface tension
  real(wp) :: geo_factor
  !! Geometric scalar factor between neighbor sectional bins.


contains

  subroutine read_file(filename)
    use tomlf, only: toml_table, toml_parse, toml_error
    
    character(*), intent(in) :: filename
    type(toml_table), allocatable :: table
    integer :: io
    type(toml_error), allocatable :: error
    open(file=filename, newunit=io, status="old")
    call toml_parse(table, io, error)
    close(io)
    if (allocated(error)) then
       print '(a)', "Error: "//error%message
       error stop
    end if

    call init(table)
  end subroutine read_file


  subroutine init(table)
    !! Reads the soot parameters from a table of information.
    !! The table is in a toml_table format, which is typically
    !! read from an input file. You can construct this table programmatically,
    !! which is good for testing.
    use tomlf, only: toml_table, toml_array, get_value, len

    type(toml_table), intent(inout) :: table
    !! A table containing the parameters
    logical :: reverse
    integer :: ival

    call get_value(table, "ND", ND)
    call get_value(table, "bins_per_2", bins_per_2)
    call get_value(table, "total_time", total_time)
    call get_value(table, "n_step", n_step)
    call get_value(table, "w_loss", w_loss)
    call get_value(table, "scg_loss", scg_loss)
    call get_value(table, "dilution", dilution)
    call get_value(table, "sat_pressure", sat_pressure)
    call get_value(table, "const_monomer", const_monomer)
    call get_value(table, "use_coag", use_coag)
    call get_value(table, "calc_coefficients", calc_coefficients)
    call get_value(table, "calc_sink_params", calc_sink_params)
    call get_value(table, "grow_bey_bound", grow_bey_bound)
    call get_value(table, "kernel", kernel)
    call get_value(table, "RR", RR)
    call get_value(table, "temperature", temperature)
    call get_value(table, "pressure", pressure)
    call get_value(table, "rho", rho)
    call get_value(table, "m_mono", m_mono)
    call get_value(table, "surface_tension", surface_tension)

    call processParams
    
  end subroutine init

  subroutine processParams 
    !! Performs some additional processing after reading in the parameters.
    !! Namely, it verifies sanity of the inputs, and it comptues some derived quantities.
    use Constants, only: PI

    integer :: i
    !! Used as an iteration counter

    NS = 40*bins_per_2
    NT = ND + NS
    if (ND .eq. 0) then
       print '(a)', "You need at least one discrete bin for monomers!"
       error stop
    else if ((ND .gt. 1) .and. (NS .gt. 0)) then
       DIS_SECT = 1
    else if  ((ND .eq. 1) .and. (NS .gt. 0)) then
       DIS_SECT = 2
    else if ((ND .gt. 0) .and. (NS .eq. 0)) then
       DIS_SECT = 3
    else
       print '(a)', "You need at least ND > 0 in the input!"
       error stop
    end if

    v_mono = m_mono/rho
    d_mono = (6*v_mono/PI)**(1._wp/3._wp)

    A = PI*d_mono**2*surface_tension/1.5_wp/1.38e-23/temperature

    geo_factor = 2**(1._wp/bins_per_2)

    call allocate_bins(NT)
    call initialize_bins(DIS_SECT, ND, NS, NT, geo_factor, v_mono)

  end subroutine processParams
end module soot_params
