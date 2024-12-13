**Particle Tunneling Through a Potential**
============

**Theory**
============
When dealing with objects on a atomic scale, like an electron, atoms or molecules, the classical physics breaks and Quantum Physics comes into play.
Quantum Tunneling is a phenomenon where the atomic size objects when thrown on a barrier with less energy than the said barrier are able to cross it.
This phenomenon is mathematically correct and arises due to the Schrondinger Wave equation's existence inside the barrier, as well as on the other side of it too.
The program tracks the changes of the wavefunction, entered by the user, for different regions and plots the wavefunction and probability density of the particle.

**Mathematical Treatment**
============
Assume a Gausian wavefunction, of the form
**exp(-a * (x - x0)^2) * cos(k * x)**
This wavefunction is defined in region 1, before the particle hits the barrier. This wavefunction is oscillatory in nature due to the lack of any potential (potential outside the barrier is 0).
![image](https://github.com/user-attachments/assets/0169c049-e9ee-4102-b288-d49469987870)
This wavefunction is only valid in region 1 and will undergo changes to when it enters the potential. 
Specifically speaking it will change to **exp(-kappa*x)**, where **kappa** is a constant defined by **kappa=√(2*m*E/hbar^2)**,whihch is an exponentially decreasing function.
In the third region wavefunction changes back to original wavefunction with slight modifications to account for the matching slop.
Since for a wavefunction to be valid it, along with its first derivate, should be continuous at all points.
The slopes at the transisition points of the wavefunctions should match.
The wavefunctions in the region 2 and 3 are multiplied by a constant to match the slopes of the two functions.

**GUI**
============
![image](https://github.com/user-attachments/assets/a178e909-0e1b-434f-a98b-ff15ca55db46)
The GUI window has a plotting area, a control panel and a results section to display the results after calculation.
The control panel can be used to modify the wavefunction, enrgy of particle, width of barrier, potiential inside the barrier.
There is a button to for plotting in the position space and a one for plotting in the fourier space or momentum space.
There are buttons to calculate the expectation values of position and momentum, and also the fourier transform of the complete wavefunction (including region 1,2 and 3 and combining it into one).
The following are the position and momentum space plots for the default wavefunction.
![image](https://github.com/user-attachments/assets/01be10ec-57f7-46a0-8086-af905b93fcd6)
![image](https://github.com/user-attachments/assets/b8835f1a-df71-47b3-bffb-5e616e6497ee)
(Note: the default is a symmetric wavefunction so it will give 0.00 for positon and momentum expectation values)
