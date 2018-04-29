# magTCSsim

This is the Magellan Telecope Control System simulator for the MagAO project. Implements a tracking telescope, and a TCP/IP server for access by instruments.

## Build
Requires mxlib.

Compile with 
```
  make -f $MXMAKEFILE magTCSsim
```


## To-Do
- Need to implement a GUI (the old one is broken), possibly by interfacing with Stellarium
- Need to implement simulation commands (probably on a different port) 
