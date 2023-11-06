permittivity1= 1.0
permittivity2= 0.8

@test Eta(3, permittivity1, permittivity2)==1.0
@test Eta(4, permittivity1, permittivity2)==0.8