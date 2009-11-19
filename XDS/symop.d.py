s = open("symop.d")
ss = s.read()
s.close()
sp= ss.split("#")

for spp in sp[1:]:
    spps =  spp.split("\n")
    s1 = "%s: ((%s, %s, '%s', '%s', '%s')," % tuple(spps[0].split())
    s2 = "'%s')," % (" ".join(spps[1:]))
    print s1 + s2
