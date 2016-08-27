#!/bin/bash

set -e

ImpactParam=0
echo "Setting ImpactParam in Toycluster parameter file to: ${ImpactParam}"
perl -pi -e 's/ImpactParam.*/ImpactParam '${ImpactParam}'/g' "ic_both_free.par"
grep -n --color=auto "ImpactParam" "ic_both_free.par"

for i in {0..9}; do
    ZeroEOrbitFrac="0.$i"

    echo "Setting ZeroEOrbitFrac in Toycluster parameter file to: ${ZeroEOrbitFrac}"
    perl -pi -e 's/ZeroEOrbitFrac.*/ZeroEOrbitFrac '${ZeroEOrbitFrac}'/g' "ic_both_free.par"
    grep -n --color=auto "ZeroEOrbitFrac" "ic_both_free.par"

    ./run.sh -o -u "ic_both_free.par"
done
# NB this is not physical
for i in {0..3}; do
    ZeroEOrbitFrac="1.$i"

    echo "Setting ZeroEOrbitFrac in Toycluster parameter file to: ${ZeroEOrbitFrac}"
    perl -pi -e 's/ZeroEOrbitFrac.*/ZeroEOrbitFrac '${ZeroEOrbitFrac}'/g' "ic_both_free.par"
    grep -n --color=auto "ZeroEOrbitFrac" "ic_both_free.par"

    ./run.sh -o -u "ic_both_free.par"
done
