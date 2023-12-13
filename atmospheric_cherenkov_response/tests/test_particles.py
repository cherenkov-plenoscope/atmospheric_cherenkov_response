import atmospheric_cherenkov_response as acr


def test_all_particles():
    particles = acr.particles._all()


def test_all_sites_valid():
    particles = acr.particles._all()
    for pk in particles:
        acr.particles.assert_valid(particles[pk])
