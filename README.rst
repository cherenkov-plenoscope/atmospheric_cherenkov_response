Atmospheric Cherenkov response
==============================
|BlackStyle| |LicenseBadge|


Cosmic particles, such as gamma-rays, can be detected with the help of the atmospheric Cherenkov-method. The cosmic particles enter earth's atmosphere and induce showers of new particles which run down to ground. Among these new particles is Cherenkov-light which can be detected on ground.
This library helps to estimate how good a given instrument can detect a certain type of cosmic particle at a certain energy (response-function).

Algorithm
---------

0) Pick a threshold photon number ``T1`` where trigger curve starts rising
(for a given type of primary)

1) Generate shower such that particle direction hits ground at 0,0;
shower direction spread over large solid angle Omega (energy-dep.)
(for charged particles)
{could also pick (0,0) at some height, but I believe for ``z``=0 the photon
scatter is smallest}

2) Divide ground in grid of spacing = mirror diameter; could e.g. without
too much trouble use up to ``M`` x ``M`` = 1000 x 1000 grid cells = 70 x 70 km^2;
grid area is ``A``, grid centered on (0,0)

3) Reset photon counter for each cell

3) For each shower, shift grid randomly in ``x``,``y`` by 1/2 mirror diameter

4) Loop over shower photons
   4.1) reject photon if angle outside FOV
   4.2) for each photon, calculate grid cell index ``ix``, ``iy``
        {easy since square grid}
   4.3) calculate distance of photon from cell center;
        keep photon if distance < ``R_Mirror``
   4.4) increment photon counter for cell
   4.5) optionally save photon in a buffer

5) Loop over grid cells
   5.1) count cells with photons > ``T1``: ``N_1``
   5.2) using trigger curve for given particle type;
        calculate trigger prob. for (real) trigger
        and randomly reject events: keep ``N_2``
        {or simply use a 2nd threshold where trigger prob=0.5}
   5.3) Increment event counters by ``N_1``, ``N_2``
        Increment error counters by ``N_1^2``, ``N_2^2``

6) For detailed simulation, optionally output photons for
   few randomly selected ``T1``-triggered cells
   (up to 10 should be fine, given that
   probably only one of 10 triggers the detailed simulation)

7) Toy effective area (x solid angle): (``N_1`` event counter/``M``^2/Nevent)*``A``*``Omega``
   error = ``sqrt(error counter)`` ...
   Somewhat better effective area: ``N_2`` event counter ...
   Final eff. area: ``N1_eff`` area x fraction of events kept in detailed sim.


Coordinate system
-----------------

::

                                    | z
                                    |                               starting pos.
                                    |                                  ___---O
                                    |                            ___---    / |
                                    |                      ___---     n  /   |
                                    |                ___---         io /     |
                                    |          ___---             ct /       |
                                    |    ___---                 re /         |
                starting altitude __|_---                     di /           |
                                    |                       y- /             |
                                    | _-------__          ar /               |
                                    |-    th    |_      im /                 |
                                    |   ni        |_  pr /                   |
                                    | ze            |  /                     |
                                    |               |/                       |
                        ____________|______________/________________________ |
                       /            |            /            /            / |
                      /            /|          //            /            /  |
                    3/            / |        / /            /            /   |
                    /            /  |      /  /            /            /    |
                   /____________/___|____/___/____________/____________/     |
                  /            /    |  /    /            /            /      |
  obs. level     /            /     |/     /    grid    /            /       |
  altitude -  -2/-  -  -  -  /  -  -X-----/  <-shift y /            /        |
               /            /      /|    /            /            /         |
              /____________/______/_____/____________/____________/          |
             /            /     -|  |  /            /            /           |
            /            /      /   | /            /            /            |
          1/            /  grid     |/            /            /             |
          /            /  shift x   /            /            /              |
         /____________/____________/____________/____________/               |
        /            /            / |          /            /                |
       /            /            /  |         /            /                 |
     0/            /            /   |        /            /                  |
     /            /            /    |       /            /                   |
    /____________/____________/____________/____________/                    |
          0            1           2|             3                          |
                                    |                                  ___---O
                                    |                            ___---
                                    |                      ___--- |
                                    |                ___---        |
                                    |          ___---               |
                                    |    ___---       azimuth       |
                  sea leavel z=0    |_---__________________________/______ x
                                    /
                                   /
                                  /
                                 /
                                /
                               /
                              /
                             /
                            /
                           /
                          / y

.. |BlackStyle| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black

.. |LicenseBadge| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
