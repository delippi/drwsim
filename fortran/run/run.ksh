#!/bin/ksh

#### METHOD 1 ##########
tasks=""
#tasks="$tasks getNatureFiles000-023_0923.ksh"
#tasks="$tasks getNatureFiles024-047_0924.ksh"
#tasks="$tasks getNatureFiles048-071_0925.ksh"
#tasks="$tasks getNatureFiles072-095_0926.ksh"
#tasks="$tasks getNatureFiles096-119_0927.ksh"
#tasks="$tasks getNatureFiles120-143_0928.ksh"
tasks="$tasks getNatureFiles144-167_0929.ksh"
#tasks="$tasks getNatureFiles168-171_0930.ksh"

tasks=""
#tasks="$tasks run000-023_0923.ksh"
tasks="$tasks run024-047_0924.ksh"
#tasks="$tasks run048-071_0925.ksh"
#tasks="$tasks run072-095_0926.ksh"
#tasks="$tasks run096-119_0927.ksh"
#tasks="$tasks run120-143_0928.ksh"
#tasks="$tasks run144-167_0929.ksh"
#tasks="$tasks run168-171_0930.ksh"

#### METHOD 2 ##########
#tasks="getNatureFiles000-023_0923.ksh                    "
#tasks="getNatureFiles024-047_0924.ksh run000-023_0923.ksh"
#tasks="getNatureFiles048-071_0925.ksh run024-047_0924.ksh"
#tasks="getNatureFiles072-095_0926.ksh run048-071_0925.ksh"
#tasks="getNatureFiles096-119_0927.ksh run072-095_0926.ksh"
#tasks="getNatureFiles120-143_0928.ksh run096-119_0927.ksh"
#tasks="getNatureFiles144-167_0929.ksh run120-143_0928.ksh"
#tasks="getNatureFiles168-171_0930.ksh run144-167_0929.ksh"
#tasks="                               run168-171_0930.ksh"

home=`pwd`

for task in $tasks; do
    pre=`echo $task | cut -c 1-3`
    if [[ $pre == "get" ]]; then; echo "sbatch tools/$task"; fi
    if [[ $pre == "run" ]]; then; echo "ksh $task"; fi
done

echo "Does this look okay? (y/n)"
read ans
if [[ $ans == "y" ]]; then

for task in $tasks; do
    pre=`echo $task | cut -c 1-3`
    if [[ $pre == "get" ]]; then; cd $home/tools; sbatch $task; fi
    if [[ $pre == "run" ]]; then; cd $home;          ksh $task; fi
done
fi

if [[ $ans == "n" ]]; then
   exit 99
fi
