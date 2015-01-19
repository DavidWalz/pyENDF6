# pyENDF6

This python module provides a minimal set of function for reading tabulated data from ENDF nuclear data files (see
https://www-nds.iaea.org/public/endf/).

#### Example
```python
import ENDF6
f = open('C006.endf')
lines = f.readlines()
sec = ENDF6.find_section(lines, MF=3, MT=3)  # total cross-section
x, y = ENDF6.read_table(sec)

figure()
plot(x, y)
xlabel('Photon energy [eV]')
ylabel('Cross-section [barn]')
show()
```
