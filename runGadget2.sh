#!/usr/bin/env bash
# @BMP /bin/bash is v3; /usr/local/bin/bash is v4.3.
# Using /usr/bin/env bash uses the correct bash at all machines :-)
# File: runGadget2.sh
# Author: Timo L. R. Halbesma <timo.halbesma@student.uva.nl>
# Date created: Fri Dec 04, 2015 03:44 PM
# Last modified: Fri Apr 22, 2016 03:35 pm
#
# Description: Compile Gadget-2, run simulation, copy Gadget-2 makefile/param

# "A good strategy for doing this in practice is to make a copy
# of the whole simulation source code together with its makefile in the output
# directory of each simulation run, and then use this copy to compile the code
# and to run the simulation. The code and its settings become then a logical
# part of the output generated in the simulation directory."
# --> We do not copy the entire directory but do copy the makefile/param files.

set -e

loglevel='INFO'

# Checking for indianness should not be necessary?
test_indianness() {
    gcc checkIndianness.c -o checkIndianness 
    ind="$(./checkIndianness)"

    if [ $ind -eq 1 ]; then
        echo "This machine is little endian."
    elif [ $ind -eq 0 ]; then
        echo "This machine is big endian."
    else 
        echo "ERROR: Unknown if big endian or little endian."
        exit 1
    fi
}

_usage() {
cat <<EOF
Usage: `basename $0` <[options]> <[filename]>
Compile Gadget-2, execute FILE with Gadget-2

Options:
  -r   --restart        set Gadget-2 restart flag to continue execution
                        implies the Gadget-2 executable will not be regenerated
                        i.e. the source will not be build.
  -l   --loglevel       set loglevel to 'ERROR', 'WARNING', 'INFO', 'DEBUG'
  -h   --help           display this help and exit
  -t   --test           test indianness of machine
  -v   --version        output version information and exit

Examples:
  `basename $0` merge_clusters.c
  `basename $0` -r merge_clusters.c
EOF
}


parse_options() {
    # It is possible to use multiple arguments for a long option. 
    # Specifiy here how many are expected.
    declare -A longoptspec
    longoptspec=( [loglevel]=1 ) 

    optspec=":rlhtv-:"
    while getopts "$optspec" opt; do
    while true; do
        case "${opt}" in
            -) # OPTARG is long-option or long-option=value.
                # Single argument:   --key=value.
                if [[ "${OPTARG}" =~ .*=.* ]] 
                then
                    opt=${OPTARG/=*/}
                    OPTARG=${OPTARG#*=}
                    ((OPTIND--))    
                # Multiple arguments: --key value1 value2.
                else 
                    opt="$OPTARG"
                    OPTARG=(${@:OPTIND:$((longoptspec[$opt]))})
                fi
                ((OPTIND+=longoptspec[$opt]))
                # opt/OPTARG set, thus, we can process them as if getopts would've given us long options
                continue 
                ;;
            l|loglevel)
                loglevel=$OPTARG
                echo "The loglevel is $loglevel"
                ;;
            r|restart)
                restart_gadget=$OPTARG
                echo "Restart flag set"
                ;;
            h|help)
                _usage
                exit 2  # 2 means incorrect usage
                ;;
            t|test)
                test_indianness
                exit 0
                ;;
            v|version)
                grep "# Version:" $0
                grep "# Last modified:" $0
                echo "Copyright (C) 2015 Timo L. R. Halbesma, BSc."
                exit 0
                ;;
        esac
    break; done
    done

    # Not sure if this is needed...
    # shift $((OPTIND-1)) 
}


# 'Main'

# Uncomment if options are required
# if [ $# = 0 ]; then _usage; fi
parse_options $@

# Unpack tar archive with source code
# ( cd source ; tar xzv --strip-components=1 -f - ) < gadget-2.0.7.tar.gz 

# Copy source code and makefile to output of file
# Then compile, do some more moving around of files
example="toycluster"
cd "example_$example"
if [ ! -d Gadget2_source ]; then
    mkdir Gadget2_source
fi
cp -r ../source/Gadget2/* ./Gadget2_source/
cp "$example.Makefile" Gadget2_source/Makefile 
cd Gadget2_source
make -j8
cd ..
mv Gadget2_source/Gadget2 .

if [ ! -d out ]; then
    mkdir out
fi

# Run code
# Working with binary files: fun! Check indianness.
if [[ -z "$ind" ]]; then
    cd ..
    echo "Checking if machine is little endian or big endian:"
    test_indianness
    cd "example_$example"
fi

# Adjust parameterfile based on indianness.
if [ $ind -eq 1 ]; then
    perl -pi -e 's/bigendian.dat/littleendian.dat/g' "$example.param" 
elif [ $ind -eq 0 ]; then
    perl -pi -e 's/littleendian.dat/bigendian.dat/g' "$example.param" 
else 
    echo "ERROR: Unknown if big endian or little endian."
    exit 1
fi

# Find number of cpu's. Note this will be off by a factor 2 for hypertreading 
# cpu's, e.g. intel i7?
# ncpu=$(grep -c ^processor /proc/cpuinfo)
ncpu=$(sysctl -n hw.ncpu)  # OSX
# Actually run the code
mpiexec -np $ncpu ./Gadget2 "$example.param"

#for example in cluster galaxy gassphere lcdm_gas; do
#    echo "example_$example"/"$example.Makefile"
#done
