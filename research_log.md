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
#### Jul 19,2021


1. Check slack's jc_reprise channel. Implement metapop model in python and run dynamics to check if the whitboard claim is true. For what value of c? For rock-paper-siccors P matrix. 
2. Finish getting the equilibria of equations once they have been transformed to a zero-sum game. 
3. Check the general way of transforming any 2x2 game into a partnership/zero-sum game. 
4. Compute equilibrium, eigenvalues, and eigenvectors of the zero-sum game.5. Get Jacobian at equilibrium, and calculate eigenvalues. Check for local stability.
6. Perform simulations under the right conditions to show coexistence, but also its lack of robustness.


#### Sep 21,2021

Sometimes I obtain coexistance of more than 2 soil or plants when sampling A at random. is this a mistake of implementation, or that the coexistence conditions are being met by chance? Check this. 

#### Sep 27,2021

I have implemented Zachs feedback about sampling real numbers instead of integers when sampling A and B without constraints, in order to ensure non-coexistence. Now its a lot better, I only get a few more-than-two plant communities. What I am going to do to fix this is: 

1. Integrate dynamics
2. Check if soil goes extinct before plant. If that happens, discard simulation
3. Check if I have reached equilibrium
2. If not, re-integrate, and get rid of extinct species (with a tolerance of tol = 1e-9)
3. Reintegrate dynamics without those species until there is only 2 species.

Non-parametrized case working correctly after implementing algorithm that ensures non-singular and with feasible equilibrium A matrices. 
Parametrized case yielding either full coexistence, full extinction or (less frequently) something in between. To avoid this, try implementing the procedure dtailed above in steps 1-3.


#### Sep 28,2021

Todos for the afternoon
1. Collect n_s simulations with a while instead of a for loop
2. Fix the detection of repeated last elements to allow for a tolerance (tol_peak). 


I am almost done implementing the equilibrium check.

Things todo for tomorrow are: 

1. Figure out what the fuck is going on with divergent dinamics all of a suddeni its because it diverges on the second integration DONE
2. Set a maximum number of integrations DONE
3. Get rid of n_s, n_p, since its always the same for me. DONE

#### Sep 29,2021

Things todo for tomorrow

1. Find peaks for plant and soil solutions DONE
2. Only if peaks of all solutions are constant, step to the next if
3. The next if is going to check whether the time between peaks is roughly constant

Actually, I think I can implement 3 when I check 2, and then do a big iff that forces all peaks of all solutions to remain equally separated and at a constant height. 

#### Sep 30,2021

Need to finish implementing cluster algorithm

Stefano's talk:

He told me that I can check for equilibrium if I evaluate the function H'(x, y) in page 130 of H-S book, and it is 0. For that, I need to find a c that satisfies this. I don't know if I can do that, because in the book they say that this constant of motion is only constant if 11.5 is valid for some c in R. 

Todos thing for tomorrow

1. Find peaks for all abundance curves
2. When checking if gaps are 0, check that there are more than just 1 (maybe a percentage of the total number of peaks)
3. Think a way of implementing Stefano's method. 

#### Oct 04,2021

Turns out that I can actually check what Stefano suggested by checking H', because the system colapses at equilibrium to a rescaled zero sum game (for c<0) or partnership game (for c>0). In the first case, I could find a negative constant such that H' = 0, and in the other case, the consant would be positive. In the second case I would not have 2-species coexistence because I would be in a saddle point, so it is unstable. I could check that I am in this point because with different initial conditions I get different outcomes (priority effects). The other reason why I would get only one species is because the matrix A is unfeasible, but in my code I am already selecting for A matrices that are feasible so I should not see deterministic dynamics for the 1-sp case when I vary initial conditions.

#### Oct 06,2021

Stefano has proved an easier way to transform one game into another by assuming that one matrix (B) is diagonal. So I need to tell Zach. This is without loss of generality, because then I can always redefine soil and/or plants to be a composition of the different soils/ plants such that one of the matrix becomes symetric.

Also, for detecting equilibria: I can integrate until I have 2 species, and then solve my system of equations to see if I am in the case of partnership or zero sum. If I am in zero sum, that is an equilibrium, and if I am in partnership, then I keep integrating to see that one species goes extinct.

#### Oct 18,2021

Draw parameters
Integrate for n timesteps
check number of species
if more than 2, keep integrating
if 2 or less, stop. 
Record number of cycles.
If after n integration cycles there is no convergence, raise an error and save
the parameters of the system (that is, matrix A and B)

#### Feb 08,2022

Check for feasibility accordingto what specified in line 144.
