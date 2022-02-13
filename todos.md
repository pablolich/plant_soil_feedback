

#### Sep 27,2021

Write down the explanation of why the matrices that abide Hofbauer and Sigmund condition are always feasible, using the Sherman-Morrison formula to demonstarte theat the sign depends on ly on mu, which is chosen always to be positive.

#### Oct 04,2021

Turns out that I can actually check what Stefano suggested by checking H', because the system colapses at equilibrium to a rescaled zero sum game (for c<0) or partnership game (for c>0). In the first case, I could find a negative constant such that H' = 0, and in the other case, the consant would be positive. In the second case I would not have 2-species coexistence because I would be in a saddle point, so it is unstable (Exercise 11.2.3, for c>0, 11.8 is a lyapunov function for 11.2). I could check that I am in this point because with different initial conditions I get different outcomes (priority effects). The other reason why I would get only one species is because the matrix A is unfeasible, but in my code I am already selecting for A matrices that are feasible so I should not see deterministic dynamics for the 1-sp case when I vary initial conditions.

Think of a way to solve Exercise 11.2.14 (page 132) of HS book. In a simpler way than transforming both A and B such that A = B^T. Maybe through only transforming B or A!

#### Oct 18,2021

I am not checking if when the system reaches 2 species, it can keeep going down to 1, or stays at 2, with the 'smart trick'. Need to implement that

#### Oct 21,2021

Erase unnecessary packages from `plant_feedback_model.py`

#### Oct 25,2021

modify integrate_PSF so that it performs multiple cycles of integration until
convergence, or until max_cycle is reached. 
Can make it recursive!
Change function of equilibrium by not making it depend on the number of species left, but the number of species whose abundances are 0!

#### Feb 09,2022

Integrate for several epsilons on each simulation that ends in 1 or 2 species while checking for convergence. Add one more loop layer. Only get simulations where three epsilon converge. figure out what to do with soil ploting. Plot with transparency/dashed
