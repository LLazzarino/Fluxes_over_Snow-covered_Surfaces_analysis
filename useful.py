import pandas as pd
import numpy as np
import os
    
def read_day(filename):

    # Impostazioni iniziali
    zs = [1, 2, 5, 10, 15, 20, 30]  # livelli FLOSS
    uidx = np.arange(0, 28, 4)  # indici u
    vidx = uidx + 1  # indici v
    widx = uidx + 2  # indici w
    Tidx = uidx + 3  # indici T

    # Percorso del file per un singolo giorno (ad esempio, apr01.dat)
    file_name = filename
    directory = ""
    file_path = os.path.join(directory, file_name)

    # Estrazione del giorno e del mese dal nome del file
    ggiorno = int(file_name[3:5])
    gmese = file_name[0:3]
    mesi = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
    gmese = mesi.index(gmese) + 1  # Trovare l'indice del mese (1-12)

    # Lettura del file
    with open(file_path, 'r') as file_read:
        # Sostituzione di '*******' con 'NA' (valori mancanti)
        file_data = file_read.read().replace('*******', ' NA')

    # Scrittura temporanea dei dati modificati in un file temporaneo
    temp_file = "temp_file.txt"
    with open(temp_file, 'w') as temp:
        temp.write(file_data)

    # Lettura dei dati come DataFrame (simile a fread in R)
    dati_giorno = pd.read_csv(temp_file, header=None, delim_whitespace=True, nrows=2160000)

    # Rinominare le colonne per u, v, w, T
    colonne = dati_giorno.columns.tolist()  # Prendere le colonne come lista

    for idx, z in zip(uidx, zs):
        colonne[idx] = f"u_{z}"
    for idx, z in zip(vidx, zs):
        colonne[idx] = f"v_{z}"
    for idx, z in zip(widx, zs):
        colonne[idx] = f"w_{z}"
    for idx, z in zip(Tidx, zs):
        colonne[idx] = f"T_{z}"

    # Assegnare la nuova lista di colonne al DataFrame
    dati_giorno.columns = colonne

    # Convertire la temperatura da Celsius a Kelvin
    dati_giorno.iloc[:, Tidx] += 273.16

    # Rimuovere il file temporaneo
    os.remove(temp_file)

    # Stampare per verificare (opzionale)
    print(dati_giorno.head())
    return dati_giorno

# Funzione per eseguire la rotazione singola/doppia/tripla
def Rotation_1_2_3(datiorari, zero_rot=True, first_rot=True, second_rot=True): #SECOND ROT IS NOW TRUE
    # zero_rot=TRUE    # to perform first rotation (<v>=0)
    # first_rot=TRUE   # to perform second rotation (<w>=0)
    # second_rot=TRUE  # to perform third rotation (<v'w'>=0)
    # datiorari: dataframe con colonne (u, v, w, T)
    
    # Estrai le componenti u, v, w
    uvw_old = datiorari.iloc[:, 0:3].values
    ndati = len(datiorari)
    
    # Filtra i dati validi (non NA)
    uvw_validi = uvw_old[~np.isnan(uvw_old).all(axis=1)]
    ndativalidi = len(uvw_validi)
    
    u_old = uvw_old[:, 0]
    v_old = uvw_old[:, 1]
    w_old = uvw_old[:, 2]
    
    control_NA = datiorari.iloc[:, 0:3].isna().sum()
    
    if ndativalidi > 1 and ndativalidi and control_NA.max() != len(u_old):
        # Valori medi
        umed = np.nanmean(u_old)
        vmed = np.nanmean(v_old)
        wmed = np.nanmean(w_old)
        
        def covariance(x, y):
        # Usa np.cov ignorando i valori NaN
            return np.nan_to_num(np.cov(x, y, rowvar=False))[0, 1]

        # Inizializzazione della matrice di covarianza 3x3
        cov0R = np.zeros((3, 3))

        # Calcolo delle covarianze tra le componenti u, v, w
        cov0R[0, 0] = covariance(u_old, u_old)
        cov0R[0, 1] = covariance(u_old, v_old)
        cov0R[0, 2] = covariance(u_old, w_old)

        cov0R[1, 0] = cov0R[0, 1]  # Simmetrico
        cov0R[1, 1] = covariance(v_old, v_old)
        cov0R[1, 2] = covariance(v_old, w_old)

        cov0R[2, 0] = cov0R[0, 2]  # Simmetrico
        cov0R[2, 1] = cov0R[1, 2]  # Simmetrico
        cov0R[2, 2] = covariance(w_old, w_old)
        
        # Generazione della matrice di rotazione
        rota = rotation_matrix_generation(umed, vmed, wmed, cov0R, zero_rot, first_rot, second_rot)
        
        # Valutazione del vettore ruotato
        urot = rota[0, 0] * u_old + rota[0, 1] * v_old + rota[0, 2] * w_old
        vrot = rota[1, 0] * u_old + rota[1, 1] * v_old + rota[1, 2] * w_old
        wrot = rota[2, 0] * u_old + rota[2, 1] * v_old + rota[2, 2] * w_old
        
        datiorariruot = pd.DataFrame({'urot': urot, 'vrot': vrot, 'wrot': wrot})
    else:
        datiorariruot = pd.DataFrame({
            'urot': [np.nan] * ndati,
            'vrot': [np.nan] * ndati,
            'wrot': [np.nan] * ndati
        })
    
    return datiorariruot

# Funzione per generare la matrice di rotazione
def rotation_matrix_generation(uu, vv, ww, cov0R, zero_rot=True, first_rot=True, second_rot=True):
    # Inizializzazione della matrice unitaria
    aa = np.eye(3)
    rota = aa.copy()
    
    if zero_rot or first_rot or second_rot:
        voriz = np.sqrt(uu**2 + vv**2)
        
        # Prima rotazione
        if zero_rot:
            alfa = np.arctan2(vv, uu) if uu != 0 and vv != 0 else 0.
            sinal = np.sin(alfa)
            cosal = np.cos(alfa)
            bb = np.array([
                [cosal, sinal, 0],
                [-sinal, cosal, 0],
                [0, 0, 1]
            ])
            rota = np.dot(bb, rota)
        
        # Seconda rotazione
        if first_rot:
            gamma = np.arctan2(ww, voriz) if ww != 0 and voriz > 0 else 0.
            singa = np.sin(gamma)
            cosga = np.cos(gamma)
            cc = np.array([
                [cosga, 0, singa],
                [0, 1, 0],
                [-singa, 0, cosga]
            ])
            rota = np.dot(cc, rota)
        
        # Terza rotazione
        if second_rot:
            cov2R = np.dot(np.dot(rota, cov0R), rota.T)
            vp2 = cov2R[1, 1]
            vpwp = cov2R[1, 2]
            wp2 = cov2R[2, 2]
            
            fi = 0.5 * np.arctan2(2 * vpwp, vp2 - wp2) if vpwp != 0 or (vp2 - wp2) != 0 else 0.
            cosfi = np.cos(fi)
            sinfi = np.sin(fi)
            dd = np.array([
                [1, 0, 0],
                [0, cosfi, sinfi],
                [0, -sinfi, cosfi]
            ])
            rota = np.dot(dd, rota)
    
    return rota


# Funzione per dividere un DataFrame in blocchi
def split_dataframe(df, row_chunks, col_chunks):
    rows_per_chunk = df.shape[0] // row_chunks
    cols_per_chunk = df.shape[1] // col_chunks
    df_matrix = []

    for i in range(0, df.shape[0], rows_per_chunk):
        row_list = []
        for j in range(0, df.shape[1], cols_per_chunk):
            chunk = df.iloc[i:i+rows_per_chunk, j:j+cols_per_chunk]
            row_list.append(chunk)
        df_matrix.append(row_list)

    return df_matrix


# Funzione per dividere il DataFrame con chiavi diverse per righe e colonne
def split_dataframe_dict(df, row_chunks, col_chunks, row_keys, col_keys):
    rows_per_chunk = df.shape[0] // row_chunks
    cols_per_chunk = df.shape[1] // col_chunks
    df_dict = {}

    for i in range(0, df.shape[0], rows_per_chunk):
        for j in range(0, df.shape[1], cols_per_chunk):
            # Creare chiave combinata tra righe e colonne
            key = f"{row_keys[i // rows_per_chunk]}_{col_keys[j // cols_per_chunk]}"
            df_dict[key] = df.iloc[i:i+rows_per_chunk, j:j+cols_per_chunk]

    return df_dict


#Plot u,v,w,T
import matplotlib
import matplotlib.pyplot as plt
#parameter to use Latex in matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
#a plot label should be easy to read
plt.rcParams.update({'font.size': 14})
import re
def plotHour(dataframe, hour= "h20", height = "z1", timetoplots = [0,-1]) :
    hour_loc = hour + "_" + height
    data = {
        'u'       : dataframe[hour_loc].iloc[:,0].values, 
        'v'       : dataframe[hour_loc].iloc[:,1].values, 
        'w'       : dataframe[hour_loc].iloc[:,2].values, 
        'T'       : dataframe[hour_loc].iloc[:,3].values
    }
    
    h = re.findall(r'\d+', hour)
    h = int(h[0])
    
    frequenza = 60  # 60 dati al secondo
    tempo_totale_secondi = 3600  # Ad esempio, 10 secondi
    time_index = pd.date_range(start=f'{h}:00:00', periods=frequenza * tempo_totale_secondi, freq='16.666666ms')


    quantities = ['u','v','w','T']

    ylabels = ['u [m/s]','v [m/s]','w [m/s]','T [K]']

    colors = ['darkgoldenrod','mediumseagreen','rebeccapurple', 'coral']

    fig,ax = plt.subplots(2,2,figsize=(12,6))
    axs = ax.flatten()

    #plot 
    for i in range(4):
        mean_value = data[quantities[i]].mean()  # Calculate mean
        var_value = data[quantities[i]].var()    # Calculate variance
        
        axs[i].plot(time_index.time,data[quantities[i]], color=colors[i], 
                    label=rf'{quantities[i]} $\quad$ [$\langle {quantities[i]} \rangle$={mean_value:.2f},$\,$$\sigma^2$$_{quantities[i]}$$_{quantities[i]}$={var_value:.2f}]')
        axs[i].grid(True)
        axs[i].legend()
        #axs[i].set_xbound(20:00:00,21:00:00)
        axs[i].set_xlim(time_index[timetoplots[0]].time(), time_index[timetoplots[1]].time())
        axs[i].set_xlabel("time")
        axs[i].set_ylabel(ylabels[i])
        #plt.xticks(rotation=45)
        #axs[i].set_xticks(df['Tempo'][::600])
        #axs[i].set_xticklabels(df['Tempo'][::600], rotation=45)
        axs[i].tick_params(axis='x', labelrotation=30)
        
    fig.suptitle(hour_loc)
   # fig.supxlabel('time')
    plt.subplots_adjust(left=0.08, bottom=0.1, right=0.9, top=0.93, wspace=0.3, hspace=0.48)
    plt.show()
    

def rotateHour(df, hour= "h5", height = "z30"):
    hour_loc = hour + "_" + height
    temp_rotated = Rotation_1_2_3(df[hour_loc])
    df[hour_loc].iloc[:, :3] = temp_rotated.values
    return df[hour_loc]