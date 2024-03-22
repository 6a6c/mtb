##Mtb Master Equation Code
### Jake Looney
### jlooney314@gmail.com or jlooney6@vols.utk.edu

This repo contains the python code for creating and solving the master equation for the 2D-birth death process modeling replication clock plasmid Mtb.

To run these files, you'll need a handful of libraries, including `scipy`, `numpy`, `cupy`, and maybe `matplotlib` and some others.

These files will need to be edited to change the input/output directories to fit your needs. Hopefully soon I'll be able to go through and add
command line arguments to make this easier.

I recommend starting with `python3 mtb-ivp.py`. This will solve the IVP on the CPU. Again, make sure you edit to make it output where and when you want it to.

If you have cuda on your machine, you can install [`cupy`](https://cupy.dev/) to run the hardware accelerated version.
I've been using cuda version 11.8, although the others should (hopefully) work.
Note that on many clusters you'll likely have to do something like `module add cuda/11.8.0-binary` or something to add cuda.
Then, you can run `python3 mtb_cuda-ivp.py` to solve the IVP much faster. Note that this is made available by [this lovely repo](https://github.com/yor-dev/cupy_ivp/tree/master), 
which makes scipy's solve\_ivp available to cupy.

This repo is still very dirty, and I'll be cleaning it up more and more soon. 
If you have any questions, feel free to email me or text me at 615-487-7162.


