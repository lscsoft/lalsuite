assert() {
    if [ "X$1" != "X$2" ]; then
        echo "ERROR: assertion failed: '$1' != '$2'"
        exit 1
    fi
}

assert \
    "$(lalpulsar_PublicSFTDirs 'H-1_H1_1800SFT_O4RUN+R1+Choft+WTKEY5-1257800000-1800.sft')" \
    'H1_1800SFT_O4RUN+R1+Choft+WTKEY5_BROADBAND-1257'
assert \
    "$(echo '/path/to/H-1_H1_1800SFT_O4RUN+R1+Choft+WTKEY5-1257800000-1800.sft' | lalpulsar_PublicSFTDirs)" \
    'H1_1800SFT_O4RUN+R1+Choft+WTKEY5_BROADBAND-1257'
assert \
    "$(lalpulsar_PublicSFTDirs -F 'mv {orig} {dir}/{name}' '/path/to/H-1_H1_1800SFT_O4RUN+R1+Choft+WTKEY5-1257800000-1800.sft')" \
    'mv /path/to/H-1_H1_1800SFT_O4RUN+R1+Choft+WTKEY5-1257800000-1800.sft H1_1800SFT_O4RUN+R1+Choft+WTKEY5_BROADBAND-1257/H-1_H1_1800SFT_O4RUN+R1+Choft+WTKEY5-1257800000-1800.sft'

assert \
    "$(lalpulsar_PublicSFTDirs 'H-5_H1_1800SFT_O5SIM+R2+Choft+WHANN_NBF0010Hz0W0008Hz0-1257800000-9000.sft')" \
    'H1_1800SFT_O5SIM+R2+Choft+WHANN_NARROWBAND-NBF0010Hz0W0008Hz0'
assert \
    "$(lalpulsar_PublicSFTDirs '/path/to/H-5_H1_1800SFT_O5SIM+R2+Choft+WHANN_NBF0010Hz0W0008Hz0-1257800000-9000.sft')" \
    'H1_1800SFT_O5SIM+R2+Choft+WHANN_NARROWBAND-NBF0010Hz0W0008Hz0'
assert \
    "$(echo '/path/to/H-5_H1_1800SFT_O5SIM+R2+Choft+WHANN_NBF0010Hz0W0008Hz0-1257800000-9000.sft' | lalpulsar_PublicSFTDirs -F 'mv {orig} {dir}/{name}')" \
    'mv /path/to/H-5_H1_1800SFT_O5SIM+R2+Choft+WHANN_NBF0010Hz0W0008Hz0-1257800000-9000.sft H1_1800SFT_O5SIM+R2+Choft+WHANN_NARROWBAND-NBF0010Hz0W0008Hz0/H-5_H1_1800SFT_O5SIM+R2+Choft+WHANN_NBF0010Hz0W0008Hz0-1257800000-9000.sft'
