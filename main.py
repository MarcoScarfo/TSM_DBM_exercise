import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import os

# LETTURA FILE E CREAZIONE STRUTTURA DATI
file = '\Data\Inconel 718\Inc718-Cast_001 - 200 s-1 - 20C\Cast-200-20C-001-in-out-signals.txt'
# file = '\Data\Inconel 718\Inc718-Cast_007 - 800 s-1 - 20C\Cast-800-20C-007-in-out-signals.txt'
# file = '\Data\Inconel 718\Inc718-Cast_050 - Statica - 20C\Inc718-Cast_50.txt'

file_path = os.path.dirname(os.path.realpath(__file__)) + file
data_frame = pd.read_csv(file_path, sep="\t")

# PARAMETRI DEL CAMPIONE
d_provino = 3  # Diametro provino
A_provino = (d_provino/2)**2*np.pi  # Area provino
E_provino = 2e5  # Modulo elastico

# STRESS
data_frame['Stress Input'] = data_frame['Input Load [N]']/A_provino
data_frame['Stress Output'] = data_frame['Output Load [N]']/A_provino

# STRAIN
data_frame['Strain Input'] = data_frame['Stress Input']/E_provino
data_frame['Strain Output'] = data_frame['Stress Output']/E_provino
data_frame['Strain Input %'] = data_frame['Strain Input']*100
data_frame['Strain Output %'] = data_frame['Strain Output']*100

# GRAFICO COMPLETO COME IMMAGINE
plt.Figure()
plt.plot(data_frame['Time [s]'], data_frame['Stress Input'])
plt.plot(data_frame['Time [s]'], data_frame['Stress Output'])
plt.show()

# FILTRAGGIO (PUNTO DI INIZIO E DI FINE)
tempo = data_frame['Time [s]'].values
stress_in = data_frame['Stress Input'].values
stress_out = data_frame['Stress Output'].values

# GRAFICO COMPLETO INTERATTIVO NEL BROWSER
fig = px.line(x=tempo, y=stress_in, title='Life expectancy in Canada')
fig.show()

passo = 5  # DA CAMBIARE PER OGNI CASO
limite = 200  # DA CAMBIARE PER OGNI CASO

differenze = abs(stress_out[passo:] - stress_out[:-passo])
posizioni = np.where(differenze > limite)[0]
id_inizio = posizioni[0] # OPPURE MANUALE
differenza_posizioni = np.diff(posizioni)
id_fine = differenza_posizioni[np.where(differenza_posizioni > 1)[0]][0] + id_inizio  # OPPURE MANUALE

# RIDUZIONE DEL GRAFICO
tempo = tempo[id_inizio:id_fine]
stress_in = stress_in[id_inizio:id_fine]
stress_out = stress_out[id_inizio:id_fine]

# GRAFICO RIDOTTO COME IMMAGINE
plt.Figure()
plt.plot(tempo, stress_in)
plt.plot(tempo, stress_out)
plt.show()

# GRAFICO RIDOTTO INTERATTIVO NEL BROWSER
fig = px.line(x=tempo, y=stress_in, title='Life expectancy in Canada')
fig.show()


a = 1