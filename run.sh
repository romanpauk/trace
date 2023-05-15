#!/bin/bash

CPUOFF=24,25,26,27,28,29,30,31
CPUSPEC=8-15,24-31
CPUSPEC_IRQ_OFF=0-7,16-31
CPUSPEC_IRQ_ON=0-31

CPUSET=/home/roman/projects/cpuset/cset

if [ ! -d /dev/cpuset ]; then
    echo "error: /dev/cpuset not found, execute 'mount -t cpuset none /dev/cpuset'"
    exit 1
fi

function cpu_online() {
    for i in ${CPUOFF//,/ }; do
        /bin/echo $1 > /sys/devices/system/cpu/cpu$i/online
    done
}

function cpu_boost() {
    INTEL_PSTATE=/sys/devices/system/cpu/intel_pstate/no_turbo
    CPUFREQ_BOOST=/sys/devices/system/cpu/cpufreq/boost

    if [ -f $INTEL_PSTATE ]; then
        if [ $1 -eq "1" ]; then
            /bin/echo 0 > $INTEL_PSTATE
        else
            /bin/echo 1 > $INTEL_PSTATE
        fi
    fi
    
    if [ -f $CPUFREQ_BOOST ]; then
        /bin/echo $1 > $CPUFREQ_BOOST
    fi
}

function cpu_irq() {
    for i in $(ls /proc/irq); do
        if [ -d /proc/irq/$i ]; then
            /bin/echo "$1" > /proc/irq/$i/smp_affinity_list || echo "setting IRQ $i failed"
        fi
    done
}

function cleanup {
    $CPUSET shield --reset
    cpu_irq "$CPUSPEC_IRQ_ON"
    cpu_online 1
    cpu_boost 1
    cpupower frequency-set -g powersave
}

trap cleanup EXIT

cpupower frequency-set -g performance
$CPUSET shield --cpu $CPUSPEC -k on
cpu_irq "$CPUSPEC_IRQ_OFF"
cpu_online 0
cpu_boost 0

$CPUSET shield --exec $@

