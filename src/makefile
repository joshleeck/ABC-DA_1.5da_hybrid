#================================================================================
# Compile and link options for ubuntu
#================================================================================
#CompileOpts=-O3 -fcheck=bounds
#CompileOpts=-O3
#LinkOpts=-L/usr/local/lib -lnetcdff -lfftpack5.1d
#LinkOptsEig=-L/usr/local/lib -lnetcdff -lfftpack5.1d -llapack -ltmglib -lrefblas

#================================================================================
# Compile and link options for Met Dept RACC cluster
#================================================================================
CompileOpts=-O3
LinkOpts=-L/usr/local/lib -L/storage/silver/metstudent/phd/bc012207/ABC_libs -lnetcdff -lfftpack5.1d
LinkOptsEig=-L/usr/local/lib -L/storage/silver/metstudent/phd/bc012207/ABC_libs -lnetcdff -lfftpack5.1d -llapack -ltmglib -lblas

#===============================================================================
# Compile and link options for MSS cluster
#===================================================================================
#CompileOpts=-O3
#LinkOpts=-L/scratch/leeck/ABC_libs -lnetcdff -lfftpack5.1d
#LinkOptsEig=-L/scratch/leeck/ABC_libs -lnetcdff -lfftpack5.1d -llapack -ltmglib -lblas

# Linking stage
Master_PrepareABC_InitState.out: Master_PrepareABC_InitState.o DefConsTypes.o SetOptions.o Initialise_um_data.o Initialise_dims.o Initialise_model_vars.o Read_um_data_2d.o Set_grid.o Set_ht_dep_cons.o Process_um_data.o Write_state_2d.o BoundaryMod.o Boundaries.o Functions.o Energy.o Calc_geost.o Calc_hydro.o Initialise_CVT.o Initialise_ObsSpec.o Effective_buoyancy.o Calc_vert_mom_source.o Calc_horiz_div.o Calc_horiz_vort.o Deallocate_dims.o Deallocate_um_data.o Deallocate_model_vars.o
	gfortran -o $@ $^ $(LinkOpts)

Master_RunNLModel.out: Master_RunNLModel.o DefConsTypes.o SetOptions.o Initialise_dims.o Initialise_model_vars.o Read_state_2d.o Set_ht_dep_cons.o ABC_NL_ModelDriver.o Calc_vert_mom_source.o Calc_horiz_div.o Calc_horiz_vort.o Effective_buoyancy.o ABC_NL_model.o Initialise_Averages.o Lscales_from_fft.o Magnitude_rms.o Boundaries.o Calc_geost.o Calc_hydro.o Energy.o Write_state_2d.o Functions.o Initialise_CVT.o Initialise_ObsSpec.o Deallocate_dims.o Deallocate_model_vars.o Deallocate_Averages.o
	gfortran -o $@ $^ $(LinkOpts)

Master_Linear_Analysis.out: Master_Linear_Analysis.o DefConsTypes.o SetOptions.o Initialise_CVT.o Initialise_ObsSpec.o
	gfortran -o $@ $^ $(LinkOptsEig)

Master_TestSuite.out: Master_TestSuite.o DefConsTypes.o SetOptions.o Read_state_2d.o Set_ht_dep_cons.o Initialise_model_vars.o Functions.o Boundaries.o Boundaries_adj.o Initialise_CVs.o Boundaries_CV.o Boundaries_CV_adj.o LinearBal_r.o LinearBal_r_adj.o Anbalw.o Anbalw_adj.o Helmholtz.o Helmholtz_adj.o HydroBal_b.o HydroBal_b_adj.o U_p.o U_p_adj.o U_v.o U_v_adj.o U_h.o U_h_adj.o U_trans.o U_trans_adj.o fft_real2spec.o fft_spec2real.o U_stddev.o U_p_inv.o Write_state_2d.o Helmholtz_inv.o Read_Covs.o InverseSymMat.o Write_one_field.o Initialise_CVT.o Initialise_ObsSpec.o Read_Obs.o ModelObservations_linear.o  ModelObservations_adj.o Initialise_dims.o ABC_NL_ModelDriver_DA.o ABC_NL_model.o Initialise_Averages.o ModelObservations.o Write_Obs.o DeAllocate_Obs.o Deallocate_Averages.o Deallocate_model_vars.o Deallocate_CVs.o Deallocate_dims.o Deallocate_CVT.o U_trans_alpha.o U_trans_alpha_adj.o U_h_alpha.o U_v_alpha.o U_h_alpha_adj.o U_v_alpha_adj.o SetHorizLoc_alpha.o SetVertLoc_alpha.o Fn_GC_corr.o
	gfortran -o $@ $^ $(LinkOptsEig)

Master_Calibration.out: Master_Calibration.o DefConsTypes.o SetOptions.o Initialise_um_data.o Initialise_dims.o Initialise_model_vars.o Read_um_data_2d.o Set_grid.o Set_ht_dep_cons.o Process_um_data.o Boundaries.o Write_state_2d.o ABC_NL_model.o Functions.o BoundaryMod.o Initialise_Averages.o Read_state_2d.o Helmholtz_inv.o LinearBal_r.o Calc_vert_cov1.o InverseSymMat.o Boundaries_CV.o Add_model_vars.o Div_model_cons.o Write_CV.o Subtract_model_vars.o U_p_inv.o Anbalw.o HydroBal_b.o Write_one_field.o Read_Covs.o Add_CVs.o Div_CV_cons.o Subtract_CVs.o CVT_Calibration_stddevs.o CVT_Calibration_vertcovs.o VertEigens.o CVT_Calibration_horizcovs.o Write_Covs.o Read_CV.o U_stddev_inv.o U_h_inv.o U_v_inv.o Initialise_CVT.o Initialise_CVs.o  fft_real2spec.o Initialise_ObsSpec.o Smooth.o Ensure_correlation_matrix.o Deallocate_dims.o Deallocate_um_data.o Deallocate_model_vars.o Deallocate_CVs.o Deallocate_Averages.o Deallocate_CVT.o SetCustomHorizSpec.o
	gfortran -o $@ $^ $(LinkOptsEig)

Master_MakeBgObs.out: Master_MakeBgObs.o DefConsTypes.o SetOptions.o Initialise_CVT.o Initialise_ObsSpec.o DeAllocate_Obs.o Initialise_model_vars.o Initialise_dims.o Read_state_2d.o Set_ht_dep_cons.o ABC_NL_ModelDriver_DA.o ModelObservations.o Functions.o Boundaries.o ABC_NL_model.o Initialise_Averages.o Write_state_2d.o Read_Covs.o Initialise_CVs.o U_trans.o Add_model_vars.o U_h.o U_v.o U_p.o U_stddev.o fft_spec2real.o Helmholtz.o LinearBal_r.o HydroBal_b.o Anbalw.o Boundaries_CV.o Write_Obs.o Calc_hydro.o Calc_geost.o Calc_vert_mom_source.o Calc_horiz_div.o Calc_horiz_vort.o Effective_buoyancy.o Deallocate_dims.o Deallocate_um_data.o Deallocate_model_vars.o Deallocate_CVs.o Deallocate_CVT.o Deallocate_Averages.o
	gfortran -o $@ $^ $(LinkOpts)

Master_ImpliedCov.out: Master_ImpliedCov.o DefConsTypes.o SetOptions.o Read_state_2d.o Initialise_CVs.o U_trans_adj.o U_trans_alpha_adj.o U_trans_alpha.o U_trans.o Write_state_2d.o U_v_adj.o U_h_adj.o U_p_adj.o U_h_alpha_adj.o U_v_alpha_adj.o U_stddev.o Boundaries_adj.o Anbalw_adj.o HydroBal_b_adj.o LinearBal_r_adj.o Helmholtz_adj.o U_v.o U_h.o U_p.o U_v_alpha.o U_h_alpha.o Boundaries.o Anbalw.o HydroBal_b.o LinearBal_r.o Helmholtz.o fft_real2spec.o fft_spec2real.o Initialise_CVT.o Initialise_aCVT.o Initialise_ObsSpec.o Functions.o Boundaries_CV.o Mul_CV_cons.o Set_ht_dep_cons.o Read_Covs.o Initialise_model_vars.o Add_model_vars.o Mul_model_vars.o Mul_model_cons.o Div_model_cons.o Write_CV.o Boundaries_CV_adj.o SetHorizLoc_alpha.o SetVertLoc_alpha.o Fn_GC_corr.o Deallocate_dims.o Deallocate_model_vars.o Deallocate_CVs.o Deallocate_CVT.o Deallocate_aCVT.o Initialise_dims.o
	gfortran -o $@ $^ $(LinkOptsEig)

Master_ImpliedTests.out: Master_ImpliedTests.o DefConsTypes.o SetOptions.o Read_state_2d.o Initialise_CVs.o U_v_adj.o U_h_adj.o Boundaries_adj.o U_v.o U_h.o Boundaries.o fft_real2spec.o fft_spec2real.o Initialise_CVT.o Initialise_ObsSpec.o Functions.o Boundaries_CV.o Set_ht_dep_cons.o Read_Covs.o Initialise_model_vars.o Write_CV.o Boundaries_CV_adj.o Deallocate_dims.o Initialise_dims.o Deallocate_model_vars.o Deallocate_CVs.o Deallocate_CVT.o
	gfortran -o $@ $^ $(LinkOpts)

Master_RawCov.out: Master_RawCov.o DefConsTypes.o SetOptions.o Read_state_2d.o Initialise_model_vars.o Set_ht_dep_cons.o Write_state_2d.o Initialise_CVT.o Initialise_ObsSpec.o Boundaries.o Deallocate_dims.o Initialise_dims.o Deallocate_model_vars.o
	gfortran -o $@ $^ $(LinkOpts)

Master_Assimilate.out: Master_Assimilate.o DefConsTypes.o SetOptions.o PenAndGrad.o Read_state_2d.o Initialise_model_vars.o Set_ht_dep_cons.o Write_state_2d.o Functions.o Boundaries.o Boundaries_adj.o Boundaries_CV.o Boundaries_CV_adj.o LinearBal_r.o LinearBal_r_adj.o Anbalw.o Anbalw_adj.o Helmholtz.o Helmholtz_adj.o HydroBal_b.o HydroBal_b_adj.o U_p.o U_p_adj.o U_v.o U_v_alpha.o U_v_adj.o U_v_alpha_adj.o U_h.o U_h_alpha.o U_h_adj.o U_h_alpha_adj.o U_trans.o U_trans_adj.o U_trans_alpha.o U_trans_alpha_adj.o fft_real2spec.o fft_spec2real.o U_stddev.o U_p_inv.o Helmholtz_inv.o Read_Covs.o InverseSymMat.o Read_Obs.o ModelObservations_linear.o  ModelObservations_adj.o DeAllocate_Obs.o Initialise_CVs.o Subtract_model_vars.o U_trans_inv.o Minus_CVs.o Add_pert_CVs.o Add_pert_model_vars.o Add_model_vars.o Mul_model_vars.o Mul_model_cons.o Div_model_cons.o Initialise_ObsSpec.o ABC_NL_ModelDriver_DA.o ModelObservations_ZeroPert.o Add_CVs.o Subtract_CVs.o Mul_CV_cons.o Write_one_field.o U_stddev_inv.o U_v_inv.o U_h_inv.o Initialise_CVT.o Initialise_aCVT.o ModelObservations.o ABC_NL_model.o Initialise_Averages.o Write_CV.o Write_Obs.o Energy.o Calc_geost.o Calc_hydro.o SetHorizLoc_alpha.o SetVertLoc_alpha.o Fn_GC_corr.o Deallocate_dims.o Deallocate_model_vars.o Deallocate_CVs.o Deallocate_Averages.o Initialise_dims.o Deallocate_CVT.o Deallocate_aCVT.o Apply_alpha_model_vars.o
	gfortran -o $@ $^ $(LinkOptsEig)

Master_GradTest.out: Master_GradTest.o DefConsTypes.o SetOptions.o PenAndGrad.o Read_state_2d.o Initialise_model_vars.o Set_ht_dep_cons.o Write_state_2d.o Functions.o Boundaries.o Boundaries_adj.o Boundaries_CV.o Boundaries_CV_adj.o LinearBal_r.o LinearBal_r_adj.o Anbalw.o Anbalw_adj.o Helmholtz.o Helmholtz_adj.o HydroBal_b.o HydroBal_b_adj.o U_p.o U_p_adj.o U_v.o U_v_alpha.o U_v_adj.o U_v_alpha_adj.o U_h.o U_h_alpha.o U_h_adj.o U_h_alpha_adj.o U_trans.o U_trans_adj.o U_trans_alpha.o U_trans_alpha_adj.o fft_real2spec.o fft_spec2real.o U_stddev.o U_p_inv.o Helmholtz_inv.o Read_Covs.o InverseSymMat.o Read_Obs.o ModelObservations_linear.o  ModelObservations_adj.o DeAllocate_Obs.o Initialise_CVs.o Subtract_model_vars.o U_trans_inv.o Minus_CVs.o Add_pert_CVs.o Add_pert_model_vars.o Add_model_vars.o Mul_model_vars.o Mul_model_cons.o Div_model_cons.o Initialise_ObsSpec.o ABC_NL_ModelDriver_DA.o ModelObservations_ZeroPert.o Add_CVs.o Subtract_CVs.o Mul_CV_cons.o Write_one_field.o U_stddev_inv.o U_v_inv.o U_h_inv.o Initialise_CVT.o Initialise_aCVT.o ModelObservations.o ABC_NL_model.o Initialise_Averages.o Write_CV.o Write_Obs.o Energy.o Calc_geost.o Calc_hydro.o SetHorizLoc_alpha.o SetVertLoc_alpha.o Fn_GC_corr.o Deallocate_dims.o Deallocate_model_vars.o Deallocate_CVs.o Deallocate_Averages.o Initialise_dims.o Deallocate_CVT.o Deallocate_aCVT.o Apply_alpha_model_vars.o
	gfortran -o $@ $^ $(LinkOptsEig)

Master_MakeEns.out: Master_MakeEns.o DefConsTypes.o SetOptions.o Read_state_2d.o Read_Obs.o ModelObservations_ZeroPert.o DeAllocate_Obs.o Initialise_model_vars.o Set_ht_dep_cons.o Write_state_2d.o Functions.o Boundaries.o Add_model_vars.o Mul_model_cons.o Div_model_cons.o Count_Obs.o Write_Obs_hx.o ModelObservations.o ABC_NL_model.o Initialise_Averages.o Deallocate_Averages.o Energy.o Deallocate_model_vars.o Initialise_dims.o Deallocate_dims.o
	gfortran -o $@ $^ $(LinkOptsEig)

Util_WriteMetaDataCVTfile.out: Util_WriteMetaDataCVTfile.o Deallocate_dims.o
	gfortran -o $@ $^ $(LinkOpts)

Util_SpecialDiags.out: DefConsTypes.o Util_SpecialDiags.o SetOptions.o Read_Covs.o Initialise_dims.o Initialise_CVs.o Read_state_2d.o U_p_inv.o Div_CV_cons.o Add_CVs.o Write_CV.o Set_ht_dep_cons.o Initialise_ObsSpec.o Initialise_CVT.o Boundaries.o Functions.o Boundaries_CV.o Helmholtz_inv.o Write_one_field.o HydroBal_b.o Anbalw.o LinearBal_r.o fft_real2spec.o U_trans_inv.o U_stddev_inv.o U_h_inv.o U_v_inv.o Deallocate_dims.o Initialise_model_vars.o Deallocate_model_vars.o Deallocate_CVs.o Deallocate_CVT.o
	gfortran -o $@ $^ $(LinkOptsEig)

Util_ImpliedCor.out: Util_ImpliedCor.o DefConsTypes.o SetOptions.o Read_state_2d.o Initialise_CVs.o U_trans_adj.o U_trans.o Write_state_2d.o U_v_adj.o U_h_adj.o U_p_adj.o U_stddev.o Boundaries_adj.o Anbalw_adj.o HydroBal_b_adj.o LinearBal_r_adj.o Helmholtz_adj.o U_v.o U_h.o U_p.o Boundaries.o Anbalw.o HydroBal_b.o LinearBal_r.o Helmholtz.o fft_real2spec.o fft_spec2real.o Initialise_CVT.o Initialise_ObsSpec.o Functions.o Boundaries_CV.o Set_ht_dep_cons.o Read_Covs.o Initialise_model_vars.o Write_CV.o Boundaries_CV_adj.o Deallocate_dims.o Deallocate_CVs.o Deallocate_CVT.o Initialise_dims.o Deallocate_model_vars.o
	gfortran -o $@ $^ $(LinkOpts)

Util_RawCor.out: Util_RawCor.o DefConsTypes.o SetOptions.o Read_state_2d.o Initialise_model_vars.o Set_ht_dep_cons.o Write_state_2d.o Initialise_CVT.o Initialise_ObsSpec.o Boundaries.o Deallocate_dims.o Initialise_dims.o Deallocate_model_vars.o
	gfortran -o $@ $^ $(LinkOpts)

Util_TestConvScaleDens.out: Util_TestConvScaleDens.o  DefConsTypes.o SetOptions.o Initialise_dims.o Initialise_model_vars.o Read_state_2d.o Set_ht_dep_cons.o Subtract_model_vars.o Write_state_2d.o ConvScaleDens.o Write_one_field.o Deallocate_dims.o Deallocate_model_vars.o Boundaries.o Read_state_2d.o Functions.o GaussEl.o fft_real2spec.o fft_spec2real.o
	gfortran -o $@ $^ $(LinkOpts)

Util_WaveSpectra.out: Util_WaveSpectra.o DefConsTypes.o SetOptions.o Initialise_dims.o Initialise_model_vars.o Read_state_2d.o Set_ht_dep_cons.o ABC_NL_ModelDriverWave.o Calc_vert_mom_source.o Calc_horiz_div.o Calc_horiz_vort.o Effective_buoyancy.o ABC_NL_model.o Initialise_Averages.o Lscales_from_fft.o Magnitude_rms.o Boundaries.o Calc_geost.o Calc_hydro.o Energy.o Write_state_2d.o Functions.o Initialise_CVT.o Initialise_ObsSpec.o Deallocate_dims.o Deallocate_model_vars.o Deallocate_Averages.o Write_one_field.o
	gfortran -o $@ $^ $(LinkOptsEig)






# Compiling stage
Master_PrepareABC_InitState.o: Master_PrepareABC_InitState.f90 DefConsTypes.o
	gfortran -c Master_PrepareABC_InitState.f90 $(CompileOpts)

DefConsTypes.o: DefConsTypes.f90
	gfortran -c $< $(CompileOpts)

SetOptions.o: SetOptions.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Initialise_um_data.o: Initialise_um_data.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Initialise_dims.o: Initialise_dims.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Initialise_model_vars.o: Initialise_model_vars.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Read_um_data_2d.o: Read_um_data_2d.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Set_grid.o: Set_grid.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Set_ht_dep_cons.o: Set_ht_dep_cons.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Process_um_data.o: Process_um_data.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Write_state_2d.o: Write_state_2d.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Read_state_2d.o: Read_state_2d.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

BoundaryMod.o: BoundaryMod.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Boundaries.o: Boundaries.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Boundaries_adj.o: Boundaries_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Functions.o: Functions.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Energy.o: Energy.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Calc_geost.o: Calc_geost.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Calc_hydro.o: Calc_hydro.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Master_RunNLModel.o: Master_RunNLModel.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

ABC_NL_ModelDriver.o: ABC_NL_ModelDriver.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Calc_vert_mom_source.o: Calc_vert_mom_source.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Calc_horiz_div.o: Calc_horiz_div.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Calc_horiz_vort.o: Calc_horiz_vort.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Effective_buoyancy.o: Effective_buoyancy.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

ABC_NL_model.o: ABC_NL_model.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Initialise_Averages.o: Initialise_Averages.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Lscales_from_fft.o: Lscales_from_fft.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Magnitude_rms.o: Magnitude_rms.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Master_Linear_Analysis.o: Master_Linear_Analysis.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Master_TestSuite.o: Master_TestSuite.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

LinearBal_r.o: LinearBal_r.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

LinearBal_r_adj.o: LinearBal_r_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Boundaries_CV.o: Boundaries_CV.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Boundaries_CV_adj.o: Boundaries_CV_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Initialise_CVs.o: Initialise_CVs.f90 DefConsTypes.o Functions.o
	gfortran -c $< $(CompileOpts)

Anbalw.o: Anbalw.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Anbalw_adj.o: Anbalw_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Helmholtz.o: Helmholtz.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Helmholtz_adj.o: Helmholtz_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Helmholtz_inv.o: Helmholtz_inv.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

HydroBal_b.o: HydroBal_b.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

HydroBal_b_adj.o: HydroBal_b_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_p.o: U_p.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_p_adj.o: U_p_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_v.o: U_v.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_v_adj.o: U_v_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_v_alpha.o: U_v_alpha.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_v_alpha_adj.o: U_v_alpha_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_h.o: U_h.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_h_adj.o: U_h_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_h_alpha.o: U_h_alpha.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_h_alpha_adj.o: U_h_alpha_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_trans.o: U_trans.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_trans_adj.o: U_trans_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_trans_alpha.o: U_trans_alpha.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_trans_alpha_adj.o: U_trans_alpha_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Apply_alpha_model_vars.o: Apply_alpha_model_vars.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

fft_real2spec.o: fft_real2spec.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

fft_spec2real.o: fft_spec2real.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_stddev.o: U_stddev.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_p_inv.o: U_p_inv.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Master_Calibration.o: Master_Calibration.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Calc_vert_cov1.o: Calc_vert_cov1.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

InverseSymMat.o: InverseSymMat.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Write_Covs.o: Write_Covs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Add_model_vars.o: Add_model_vars.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Mul_model_vars.o: Mul_model_vars.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Div_model_cons.o: Div_model_cons.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Mul_model_cons.o: Mul_model_cons.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Write_CV.o: Write_CV.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Subtract_model_vars.o: Subtract_model_vars.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Write_one_field.o: Write_one_field.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Read_CV.o: Read_CV.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Read_Covs.o: Read_Covs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Add_CVs.o: Add_CVs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Div_CV_cons.o: Div_CV_cons.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Mul_CV_cons.o: Mul_CV_cons.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Subtract_CVs.o: Subtract_CVs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

CVT_Calibration_stddevs.o: CVT_Calibration_stddevs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

CVT_Calibration_vertcovs.o: CVT_Calibration_vertcovs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

VertEigens.o: VertEigens.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

SetHorizLoc_alpha.o: SetHorizLoc_alpha.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

SetVertLoc_alpha.o: SetVertLoc_alpha.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

SetCustomHorizSpec.o: SetCustomHorizSpec.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Fn_GC_corr.o: Fn_GC_corr.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

CVT_Calibration_horizcovs.o: CVT_Calibration_horizcovs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Initialise_CVT.o: Initialise_CVT.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Initialise_aCVT.o: Initialise_aCVT.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_stddev_inv.o: U_stddev_inv.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_h_inv.o: U_h_inv.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_v_inv.o: U_v_inv.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Master_MakeBgObs.o: Master_MakeBgObs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Initialise_ObsSpec.o: Initialise_ObsSpec.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Count_Obs.o: Count_Obs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

DeAllocate_Obs.o: DeAllocate_Obs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Deallocate_dims.o: Deallocate_dims.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Deallocate_Averages.o: Deallocate_Averages.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Deallocate_CVs.o: Deallocate_CVs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Deallocate_CVT.o: Deallocate_CVT.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Deallocate_aCVT.o: Deallocate_aCVT.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Deallocate_um_data.o: Deallocate_um_data.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Deallocate_model_vars.o: Deallocate_model_vars.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

ABC_NL_ModelDriver_DA.o: ABC_NL_ModelDriver_DA.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

ModelObservations.o: ModelObservations.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

ModelObservations_linear.o: ModelObservations_linear.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

ModelObservations_adj.o: ModelObservations_adj.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Master_ImpliedCov.o: Master_ImpliedCov.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Master_ImpliedTests.o: Master_ImpliedTests.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Master_RawCov.o: Master_RawCov.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Read_Obs.o: Read_Obs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Write_Obs.o: Write_Obs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Write_Obs_hx.o: Write_Obs_hx.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Master_Assimilate.o: Master_Assimilate.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Master_MakeEns.o: Master_MakeEns.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

PenAndGrad.o: PenAndGrad.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

U_trans_inv.o: U_trans_inv.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Minus_CVs.o: Minus_CVs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Add_pert_CVs.o: Add_pert_CVs.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Add_pert_model_vars.o: Add_pert_model_vars.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

ModelObservations_ZeroPert.o: ModelObservations_ZeroPert.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Smooth.o: Smooth.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Ensure_correlation_matrix.o: Ensure_correlation_matrix.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Util_SpecialDiags.o: Util_SpecialDiags.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Util_ImpliedCor.o: Util_ImpliedCor.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Util_RawCor.o: Util_RawCor.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Master_GradTest.o: Master_GradTest.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Util_WriteMetaDataCVTfile.o: Util_WriteMetaDataCVTfile.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Util_TestConvScaleDens.o: Util_TestConvScaleDens.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

ConvScaleDens.o: ConvScaleDens.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

GaussEl.o: GaussEl.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

Util_WaveSpectra.o: Util_WaveSpectra.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)

ABC_NL_ModelDriverWave.o: ABC_NL_ModelDriverWave.f90 DefConsTypes.o
	gfortran -c $< $(CompileOpts)



all: Master_PrepareABC_InitState.out Master_RunNLModel.out Master_Linear_Analysis.out Master_TestSuite.out Master_Calibration.out Master_MakeBgObs.out Master_ImpliedCov.out Master_ImpliedTests.out Master_RawCov.out Master_Assimilate.out Master_GradTest.out Master_MakeEns.out Util_WriteMetaDataCVTfile.out Util_SpecialDiags.out Util_ImpliedCor.out Util_RawCor.out Util_TestConvScaleDens.out Util_WaveSpectra.out

clean:
	rm -rf *.o *.out *.mod
