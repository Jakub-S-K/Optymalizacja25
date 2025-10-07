import pandas as pd

import matplotlib.pyplot as plt

# Wczytaj dane z pliku CSV (średnik jako separator, przecinek jako separator dziesiętny)
df = pd.read_csv('symulacja_lab0.csv', sep=';', decimal=',', header=None)

# Nazwij kolumny dla czytelności (4 kolumny - ostatnia jest pusta)
df.columns = ['x', 'y1', 'y2', 'empty']

# Usuń pustą kolumnę
df = df.drop('empty', axis=1)

# Rysuj wykres
plt.plot(df['x'], df['y1'], label='y1')
plt.plot(df['x'], df['y2'], label='y2')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Wykres z pliku symulacja_lab0.csv')
plt.legend()
plt.grid(True)
plt.show()