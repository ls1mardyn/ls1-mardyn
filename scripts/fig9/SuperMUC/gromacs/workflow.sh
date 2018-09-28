# Workflow for simulation of TIP4P water in GROMACS based on:
# http://wbarnett.us/tutorials/1_tip4pew_water/
# if not stated otherwise each step needs only the files created in the steps before
# grompp always also needs the topol.top file even if not explicitly stated


#create topology file
echo\
"#include \"oplsaa.ff/forcefield.itp\"\n
#include \"oplsaa.ff/tip4pew.itp\"\n
\n
[ System  ]\n
TIP4PEW in water\n

[ Molecules  ]\n
" > topol.top

# generates structure file (conf.gro) with ~16913 Molecules of TIP4P water
# also updates topology
gmx solvate -cs tip4p -o conf.gro -box 8 8 8 -p topol.top


# THIS STEP DID NOT TERMINATE ON SUPERMUC SO I SKIPPED IT
# generate run input file for energy minimization (flexible)
# needs min_flex.mdp
#gmx grompp -f min_flex.mdp -o min_flex -pp min_flex -po min_flex
# run energy minimization (flexible)
#gmx mdrun -deffnm min_flex

# generate run input file for energy minimization (rigid)
# needs min_rigid.mdp
gmx grompp -f min_rigid.mdp -o min_rigid -pp min_rigid -po min_rigid
# run energy minimization (rigid)
gmx mdrun -deffnm min_rigid

# generate run input file for temperature equilibration
# needs eql_temp.mdp
gmx grompp -f eql_temp.mdp -o eql_temp -pp eql_temp -po eql_temp -c min_rigid -t min_rigid
# run temperature equilibration (NVT)
gmx mdrun -deffnm eql_temp

# generate run input file for pressure equilibration
# needs eql_pres.mdp
gmx grompp -f eql_pres.mdp -o eql_pres -pp eql_pres -po eql_pres -c eql_temp -t eql_temp
# run pressure equilibration (NPT)
gmx mdrun -deffnm eql_pres

# stack the equilibrated scenario three times in every dimension
# results in ~456651 molecules
gmx genconf -f eql_pres.gro -o eql.gro -nbox 3 3 3

# update topology file
NewNumMols=sed -n '2{p;q}' eql.gro
sed -i "s|\(SOL *\)[0-9]*|\1 ${NewNumMols}|g" topol.top

# generate run input file for benchmark
# needs bench.mdp
gmx grompp -f bench.mdp -o bench -pp bench -po bench -c eql
# run benchmark
gmx mdrun -deffnm bench

