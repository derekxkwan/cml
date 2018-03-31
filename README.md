# [cmlsynth~] v0.1alpha for Pure Data
(released under gpl v3.0)



Based on "Coupled Map Lattices as Musical Instruments" by Janko Gravner and Kyle Johnson, this Pd external is the product of Derek Kwan working together with Kyle Johnson. Work on this external is in the very initial (testing) stages and not reflective of initial final product.


## args
* max lattice size, specified lattice size, epsilon, omega, k
    * (epsilon, omega, and k are bounded by [0,1])

## inlets
* messages, lattice size, eps, omega, k

## flags 
* -array (followed by symbol) lets the user initialize internal buffer with array rather than noise. Size of the internal buffer used will be the minimum of array size and specified lattice size.

## messages
Possible messages into the first inlet include
* eps (float)
* n (float) - specifies lattice size
* omega (float)
* k (float)
* read (symbol) - read array
* reset - reset internal buffer and phasex
