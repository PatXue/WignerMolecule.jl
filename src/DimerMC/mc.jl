struct DimerMC
    T::Float64          # Temperature
    init_T::Float64     # Initial temperature (for thermalization)
    params::WignerParams

    spins::PeriodicMatrix{SpinVector}
    ηs::PeriodicMatrix{SpinVector}

    ηks::Array{ComplexF64, 3}       # Fourier transformed ηs

    outdir::String # Output directory for spin plots
    savefreq::Int  # No. of sweeps between saving spins
end