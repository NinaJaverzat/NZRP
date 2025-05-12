#!/bin/bash

################################################################################

# This function takes in input a string.
# If a folder with name given by the input string (locally) exists, it does nothing
# If the folder does not (locally) exists, the function (locally) creates it
function build_dir() {
	if [[ ! -d $1 ]]; then
		mkdir $1 &
		wait
	fi
}

# Check the number of tasks with name = process_name
# If they are found to be == max_tasks, sleeps five seconds and then checks again
function check_and_sleep() {

  max_tasks=$1
  process_name=$2

  running_tasks=`ps -C ${process_name} --no-headers | wc -l`
  while (($running_tasks >= $max_tasks)); do
    sleep 5
    running_tasks=`ps -C ${process_name} --no-headers | wc -l`
  done

}

################################################################################

build_dir "res"
build_dir "conv_res"

cd code

make o=1
make clear

L=$1
T=$2
TOT=$3
num_cores=$4

for i in $(seq 1 $TOT); do
  ./rp.exe $L $T $i &
  sleep 1
  check_and_sleep $num_cores rp.exe
done

wait

cd ../conv

make
make clear

./convolutor.exe $L $T $TOT &

wait
