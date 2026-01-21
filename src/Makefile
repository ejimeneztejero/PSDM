FC	= mpif90
FC2	= gfortran
FFLAGS	= -O3 -fallow-argument-mismatch
LFLAGS  = $(FFLAGS)

EXE     = psdm_run
EXE2    = modelconv
EXE3    = psdm2twt
EXE4    = tomo_grids
EXE5    = psdmsum

SRC     = modules.f90 \
          00-main.f90 \
          get_was_data.f90 get_mcs_data.f90 \
          optimization.f90 \
          calculate_PSDM_RTM.f90 \
          calculate_PSDM_eikonal.f90 \
          calculate_PSDM_kirchoff.f90 \
          dynamic_window.f90 \
          solve_eikonal_FSM.f90 \
          solver_ac_forward_ram.f90 \
          solver_ac_forward_ram_tt.f90 \
          solver_ac_backward_ram.f90 \
          solver_ac_forward_disk.f90 \
          solver_ac_backward_disk.f90 \
          solver_PML_matrices.f90 \
          solver_space_derivs_clean.f90 \
          solver_PML_coeffs.f90 \
          solver_FD_coeffs.f90 \
          auxiliar.f90 \
          auxiliar_PSDM.f90 \
          interpolation_akima.f90 akima433.f90 linear_interpolation.f90 \
          filters.f90 \
          fftpack.f90

SRC2 = modelconv_main.f90
SRC3 = psdm2twt_main.f90
SRC4 = tomo_grids_main.f90
SRC5 = psdmsum_main.f90

.SUFFIXES: .f90 .o

OBJ = $(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all: $(EXE) $(EXE2) $(EXE3) $(EXE4) $(EXE5)

$(EXE): $(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ)

$(EXE2): $(SRC2)
	$(FC2) $(FFLAGS) -o $@ $(SRC2)

$(EXE3): $(SRC3)
	$(FC2) $(FFLAGS) -o $@ $(SRC3)

$(EXE4): $(SRC4)
	$(FC2) $(FFLAGS) -o $@ $(SRC4)

$(EXE5): $(SRC5)
	$(FC2) $(FFLAGS) -o $@ $(SRC5)

clean:
	rm -f $(OBJ) $(EXE) $(EXE2) $(EXE3) $(EXE4) $(EXE5) core
