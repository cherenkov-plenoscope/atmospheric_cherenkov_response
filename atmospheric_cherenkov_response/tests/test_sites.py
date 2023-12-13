import atmospheric_cherenkov_response as acr


def test_all_sites():
    sites = acr.sites._all()


def test_all_sites_valid():
    sites = acr.sites._all()
    for sk in sites:
        acr.sites.assert_valid(sites[sk])
