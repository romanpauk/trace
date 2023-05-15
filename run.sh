#!/bin/bash

CPUOFF=4,5
CPUSPEC=4-5,10-11
CPUSPEC_IRQ_OFF=1-3,6-9
CPUSPEC_IRQ_ON=1-11

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

function cpu_boost_disable() {
    if [ -f /sys/devices/system/cpu/intel_pstate/no_turbo ]; then
        /bin/echo $1 > /sys/devices/system/cpu/intel_pstate/no_turbo
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
    cpu_boost_disable 0
    cpupower frequency-set -g powersave
}

trap cleanup EXIT

cpupower frequency-set -g performance
$CPUSET shield --cpu $CPUSPEC -k on
cpu_irq "$CPUSPEC_IRQ_OFF"
cpu_online 0
cpu_boost_disable 1

$CPUSET shield --exec $@

