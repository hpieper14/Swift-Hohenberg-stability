# Swift-Hohenberg-stability
This code is written by Hannah Pieper and accompanies the paper "The Maslov 
Index, degenerate crossings, and the stability of pulse solutions to the 
Swift-Hohenberg equation" by Margaret Beck, Jonathan Jaquette and Hannah Pieper. 
** ADD LINK TO PAPER WHEN UP ** 

## Subdirectories 
[results](results): Contains the numerics presented in Section 5 of the paper. 
Run either [main.mlx](results/main.mlx) or [main.m](results/main.m). Both files
run the same code, but the live script [main.mlx](results/main.mlx) contains 
some mathematical exposition. 

[source]{source): Contains the source code that generates the results. This 
code is organized into 3 modules: 
* [Fourier Series](source/@FourierSeries): contains code for a `FourierSeries` object 
and operations.
* [Pulse Solution](source/@PulseSolution): contains code for a `PulseSolution`
object. This approximates the pulse solution to the Swift-Hohenberg equation 
using a Fourier series and Newton's method. 
* [Conjugate Points](source/@ConjugatePoints): contains code for a `ConjugatePoints` 
object. This computes the conjugate points associated to the pulse solution. 




## Copyright and License 

Copyright (C) 2023  M Beck, J Jaquette, and H Pieper.

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
[GNU General Public License](LICENSE) for more details.