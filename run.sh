#!/bin/bash -eux

script_dir=$(cd $(dirname $BASH_SOURCE); pwd)
$script_dir/oop_nodes_genetic $@
python $script_dir/plot_timeseries.py

