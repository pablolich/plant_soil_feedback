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
2. Get rid of extinct species (with a tolerance of tol = 1e-9)
3. Reintegrate dynamics without those species until there is only 2 species.

Non-parametrized case working correctly after implementing algorithm that ensures non-singular and with feasible equilibrium A matrices. 
Parametrized case yielding either full coexistence, full extinction or (less frequently) something in between. To avoid this, try implementing the procedure dtailed above in steps 1-3.
