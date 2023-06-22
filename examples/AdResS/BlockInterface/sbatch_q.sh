#!/usr/bin/zsh
: "${E_1k:=0}"
: "${E_4k:=0}"
: "${E_16k:=0}"
: "${E_64k:=0}"
: "${E_256k:=0}"
: "${E_1024k:=0}"

if [[ $E_1k -gt 0 ]]
then
  sbatch ./CI_INIT_1k_C6H12_0.cmd
  sbatch ./CI_INIT_1k_C6H12_1.cmd
  sbatch ./CI_INIT_1k_CH4_0.cmd
  sbatch ./CI_INIT_1k_CH4_1.cmd
fi

if [[ $E_4k -gt 0 ]]
then
  sbatch ./CI_INIT_4k_C6H12_0.cmd
  sbatch ./CI_INIT_4k_C6H12_1.cmd
  sbatch ./CI_INIT_4k_CH4_0.cmd
  sbatch ./CI_INIT_4k_CH4_1.cmd
fi

if [[ $E_16k -gt 0 ]]
then
  sbatch ./CI_INIT_16k_C6H12_0.cmd
  sbatch ./CI_INIT_16k_C6H12_1.cmd
  sbatch ./CI_INIT_16k_CH4_0.cmd
  sbatch ./CI_INIT_16k_CH4_1.cmd
fi

if [[ $E_64k -gt 0 ]]
then
  sbatch ./CI_INIT_64k_C6H12_0.cmd
  sbatch ./CI_INIT_64k_C6H12_1.cmd
  sbatch ./CI_INIT_64k_CH4_0.cmd
  sbatch ./CI_INIT_64k_CH4_1.cmd
fi

if [[ $E_256k -gt 0 ]]
then
  sbatch ./CI_INIT_256k_C6H12_0.cmd
  sbatch ./CI_INIT_256k_C6H12_1.cmd
  sbatch ./CI_INIT_256k_CH4_0.cmd
  sbatch ./CI_INIT_256k_CH4_1.cmd
fi

if [[ $E_1024k -gt 0 ]]
then
  sbatch ./CI_INIT_1024k_C6H12_0.cmd
  sbatch ./CI_INIT_1024k_C6H12_1.cmd
  sbatch ./CI_INIT_1024k_CH4_0.cmd
  sbatch ./CI_INIT_1024k_CH4_1.cmd
fi