# public

This repo contains code that I want to share with the public.

If you want to know what things are for and how to use them, and the readme's aren't enough, just let me know via mladen.ivkovic [at] hotmail [dot] com


So far, it contains:


- godunov:
	a 1D simple godunov solver, mainly for pedagogical purposes
	
-  python_modules:
	Some python3 modules. Just add their location to your pythonpath and you're good to go.
	
	- meshless :
		routines to compute effective surfaces a la Hopkins 2015 or Ivanova 2013. Inefficient, meant 
		to be used for easy checks and plotting, not serious computations.
	
	- physics :
		some very common physics functions and mainly various constants. Meant to be easily imported
		during interactive sessions so you don't have to copypaste constants from wikipedia all
		the time.

- riemann: 
	Simple 1D riemann solvers, including exact and approximate solvers and vacuum situations,
	mainly for pedagogical purposes
