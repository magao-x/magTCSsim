# magTCSsim

This is the Magellan Telecope Control System simulator for the MagAO project. Implements a tracking telescope, and a TCP/IP server for access by instruments.

## Build
Requires mxlib.

Compile with 
```
  make -f $MXMAKEFILE magTCSsim
```

## Simulator Commands

These are commands used to control the simulation that aren't in the normal TCS commands. 

Connect with:
```
$ telnet localhost 5811
```

Then send commands starting with `tcssim` with:
```
$ tcssim <command> <arg>
```

quit
simulating start/stop
tracking start/stop
catobj name
catra ra
catdec dec
dome open
dome close


## To-Do
- Need to implement a GUI (the old one is broken), possibly by interfacing with Stellarium
- Need to implement simulation commands (probably on a different port) 
