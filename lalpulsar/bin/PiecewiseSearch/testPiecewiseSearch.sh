if ! ${PYTHON} -c 'import sympy' >/dev/null 2>&1; then
    echo "Skipping test: requires sympy"
    exit 77
fi

# --- GW170817 search

# run code
lalpulsar_PiecewiseSearch \
    --loglevel INFO \
    --tbankcode GW170817 \
    --SFTFiles \
    --h0=0.1 --cosi=0 --psi=0 --phi0=0 \
    --detectors=H1 --noise_sqrt_Sh=1 \
    --fmin=99.5 --fmax=100

# test for expected directories
test -d GW170817_fmin-99.5_fmax-100.0_H1_BS-0_j-0
test -d Template_Count_And_Timing_For_GW170817_fmin-99.5_fmax-100.0_H1
test -f Temps_For_GW170817_fmin-99.5_fmax-100.0_H1_i-0_tqdm_Output.txt

# --- GW190425 search

# run code
lalpulsar_PiecewiseSearch \
    --loglevel INFO \
    --tbankcode GW190425 \
    --knotnum 2 \
    --SFTFiles \
    --h0=0.1 --cosi=0 --psi=0 --phi0=0 \
    --detectors=H1 --noise_sqrt_Sh=1 \
    --fmin=99.5 --fmax=100

# test for expected directories
test -d GW190425_fmin-99.5_fmax-100.0_knotnum-2_H1_BS-0_j-0
test -d Template_Count_And_Timing_For_GW190425_fmin-99.5_fmax-100.0_knotnum-2_H1
test -f Temps_For_GW190425_fmin-99.5_fmax-100.0_knotnum-2_H1_i-0_tqdm_Output.txt
