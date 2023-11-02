#! /bin/bash
set -u
max_cmd=5
cmd_file="cmdpal.dat"
cur=0
line=1
for ((;;)) do
    if [[ "$cmd_file" == "-" ]]; then
        if ! read cmd; then
            break
        fi
    else
        if [[ $line -gt $(cat $cmd_file | wc -l) ]]; then
            break
        fi
        cmd=$(sed -n ${line}p $cmd_file)
    fi
    line=$[line+1]
    echo "$cmd"
    sh -c "$cmd" &
    pid[$cur]=$!
    cur=$[cur+1]
    while [[ $cur -ge $max_cmd ]]; do
        pid_list=""
        for i in `seq 0 1 $[max_cmd-1]`; do
            if [[ -z "$pid_list" ]]; then
                pid_list=${pid[i]}
            else
                pid_list=$pid_list,${pid[i]}
            fi
        done
        cur=0
        for i in $(ps -p $pid_list -o pid | sed 1d); do
            pid[$cur]=$i
            cur=$[cur+1]
        done
        if [[ $cur -ge $max_cmd ]]; then
            sleep 1
        fi
    done
done
wait