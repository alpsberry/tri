tri
===
A simple finite element method library<br />
Read mesh information generated by Triangle [http://www.cs.cmu.edu/afs/cs.cmu.edu/project/quake/public/www/triangle.html](http://www.cs.cmu.edu/afs/cs.cmu.edu/project/quake/public/www/triangle.html)<br /> 

Usage
--------
Read parameters from default input file

	./tri

Read parameters from myinput.input

	./tri myinput.input

Generate square mesh with normalized triangulation and refinement: first go to /meshgen, then compile in the normal way, run with the following where N is the number of columns and rows. Now it can only refine the leftmost column.

	./meshgen N

Changelog
--------
> Feb 21, 2015
* Separated the sparse linear system solver with fem solving system
> 
> Jan 25, 2015
* Rewrote definitions of most classes
* Temporarily removed the continuous FEM part which will be added after modification
> 
> May 27
* throw exceptions to gain better control of runtime error
* modified the output format in *.output
> 
> May 21
* correct the edge integral part, now deal with mesh with hanging nodes correctly
* code style formatted
>
> May 16
* DG solver is able to deal with dirichlet boundary condition
>
> May 15
* add the part of computing error in H1 norm
* remove template, clean some code and add some comments
>
> May 14
* fix the part of computing error in L2 norm
>
> Apr. 28
* the Discontinous Galerkin solving system can now deal with refinement with hanging nodes
* a simple tool to generate normalized mesh and refinement is provided, see usage
* users could specify refinement times and need to provide refinement information as is generated by the tool, meshgen, in *.ref0
>
> Apr. 21
* added the function of the exact solution of the problem
* added the part of computing error
* added the parameter to control whether to compute and output error or not
>
> Apr. 19
* modified all files with c++ template in order to deal with different problem definition for different schemes
* added a matlab visulization code plotTri.m
* added some sample mesh info files: dat/square.1.node, dat/square.1.ele, dat/square.1.edge
>
> Apr. 18
* Added a Discontinuous Galerkin solving system
* a sample parameter file for DG tri.input, and for FEM triFEM.input
>
> Mar. 28
* seperate BasicSolvingSystem and MySolvingSystem, rename the latter FEMSolvingSystem
>
> Mar. 27
* add some parameters to gain control of input/output
* a sample parameter setting file is tri.input
* able to choose between UMFPACK and SuperLU to solve the spare matrix
>
> Mar. 26
* the UMFPACK solving part is now working
* added SuperLU solving part
* rearranged some structure
* added some comments
>
> Mar. 25
* version 0.01
* known problem: UMFPACK solving part does not work
