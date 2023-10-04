##############################
Atmospheric Cherenkov response
##############################

|TestStatus| |PyPiStatus| |BlackStyle| |PackStyleBlack| |LicenseBadge|


Cosmic particles, such as gamma-rays, can be detected with the help of the atmospheric Cherenkov-method. The cosmic particles enter earth's atmosphere and induce showers of new particles which run down to ground. Among these new particles is Cherenkov-light which can be detected on ground.
This library helps to estimate how good a given instrument can detect a certain type of cosmic particle at a certain energy (response-function).


*********
Algorithm
*********


1.  Pick a threshold photon number ``T1`` where trigger curve starts rising
    (for a given type of primary)

2.  Generate shower such that particle direction hits ground at 0,0;
    shower direction spread over large solid angle Omega (energy-dep.)
    (for charged particles)
    {could also pick (0,0) at some height, but I believe for ``z`` =0 the photon
    scatter is smallest}

3.  Divide ground in grid of spacing = mirror diameter; could e.g. without
    too much trouble use up to ``M`` x ``M`` = 1000 x 1000 grid cells = 70 x 70 km^2;
    grid area is ``A``, grid centered on (0,0)

4.  Reset photon counter for each cell

5.  For each shower, shift grid randomly in ``x``, ``y`` by 1/2 mirror diameter

6.  Loop over shower photons

    a.  reject photon if angle outside FOV
    b.  for each photon, calculate grid cell index ``ix``, ``iy``
        (easy since square grid)
    c.  calculate distance of photon from cell center;
        keep photon if distance < ``R_Mirror``
    d.  increment photon counter for cell
    e.  optionally save photon in a buffer

7.  Loop over grid cells

    a.  count cells with photons > ``T1`` : ``N_1``
    b.  using trigger curve for given particle type;
        calculate trigger prob. for (real) trigger
        and randomly reject events: keep ``N_2``
        (or simply use a 2nd threshold where trigger prob=0.5)
    c.  Increment event counters by ``N_1`` , ``N_2``
        Increment error counters by ``N_1`` ^2, ``N_2`` ^2

8.  For detailed simulation, optionally output photons for
    few randomly selected ``T1`` -triggered cells
    (up to 10 should be fine, given that
    probably only one of 10 triggers the detailed simulation)

9.  Toy effective area (x solid angle): (``N_1`` event counter/``M`` ^2 / Nevent)* ``A`` * ``Omega``
    error = ``sqrt(error counter)`` ...
    Somewhat better effective area: ``N_2`` event counter ...
    Final eff. area: ``N1_eff`` area x fraction of events kept in detailed sim.


*****************
Coordinate system
*****************


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

.. |TestStatus| image:: https://github.com/cherenkov-plenoscope/atmospheric_cherenkov_response/actions/workflows/test.yml/badge.svg?branch=main
    :target: https://github.com/cherenkov-plenoscope/atmospheric_cherenkov_response/actions/workflows/test.yml

.. |PyPiStatus| image:: https://img.shields.io/pypi/v/atmospheric_cherenkov_response_cherenkov-plenoscope-project
    :target: https://pypi.org/project/atmospheric_cherenkov_response_cherenkov-plenoscope-project

.. |PackStyleBlack| image:: https://img.shields.io/badge/pack%20style-black-000000.svg
    :target: https://github.com/cherenkov-plenoscope/black_pack

.. |LicenseBadge| image:: https://img.shields.io/badge/License-MIT-yellow.svg
    :target: https://opensource.org/licenses/MIT
