def test_biofuel_potential_not_exceeded(carrier_prod, biofuel_potential, location):
    biofuel_generation = carrier_prod.sel(techs="biofuel", locs=location).sum("timesteps").item()
    assert biofuel_generation <= biofuel_potential
