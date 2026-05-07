# DEFINIZIONE DELLE SIMULAZIONI DA ESEGUIRE
Simulations = [
    #(Re = 4600, max_order = 3, Lmm = "[201, 263]", Lmmcj = "[560, 622]"),
    #(Re = 4600, max_order = 5, Lmm = "[201, 263]", Lmmcj = "[560, 622]"),
    #(Re = 4600, max_order = 7, Lmm = "[201, 263]", Lmmcj = "[560, 622]"),
    
    #(Re = 4700, max_order = 3, Lmm = "[201, 269]", Lmmcj = "[565, 630]"),
    #(Re = 4700, max_order = 5, Lmm = "[201, 269]", Lmmcj = "[565, 630]"),
    #(Re = 4700, max_order = 7, Lmm = "[201, 269]", Lmmcj = "[565, 630]"),

    #(Re = 4800, max_order = 3, Lmm = "[201, 271]", Lmmcj = "[562, 632]"),
    #(Re = 4800, max_order = 5, Lmm = "[201, 271]", Lmmcj = "[562, 632]"),
    #(Re = 4800, max_order = 7, Lmm = "[201, 271]", Lmmcj = "[562, 632]"),
    
   # (Re = 4900, max_order = 3, Lmm = "[201, 276]", Lmmcj = "[565, 639]"),
    #(Re = 4900, max_order = 5, Lmm = "[201, 276]", Lmmcj = "[565, 639]"),
    #(Re = 4900, max_order = 7, Lmm = "[201, 276]", Lmmcj = "[565, 639]"),

    #(Re = 5000, max_order = 3, Lmm = "[201, 278]", Lmmcj = "[565, 642]"),
    #(Re = 5000, max_order = 5, Lmm = "[201, 278]", Lmmcj = "[565, 642]"),
    #(Re = 5000, max_order = 7, Lmm = "[201, 278]", Lmmcj = "[565, 642]"),
    
    #(Re = 5100, max_order = 3, Lmm = "[201, 281]", Lmmcj = "[567, 647]"),
    #(Re = 5100, max_order = 5, Lmm = "[201, 281]", Lmmcj = "[567, 647]"),
    #(Re = 5100, max_order = 7, Lmm = "[201, 281]", Lmmcj = "[567, 647]"),

    (Re = 5200, max_order = 3, Lmm = "[201, 288]", Lmmcj = "[567, 654]"),
    (Re = 5200, max_order = 5, Lmm = "[201, 288]", Lmmcj = "[567, 654]"),
    (Re = 5200, max_order = 7, Lmm = "[201, 288]", Lmmcj = "[567, 654]"),
    
]

path_input = "input/U_channel.jl"

# Crea la cartella principale per i risultati se non esiste
mkpath("results")

# =====================================================================
# CREAZIONE DEL LAUNCHER TEMPORANEO
# =====================================================================
launcher_code = """
using SparseArrays, LinearAlgebra, MAT, DelimitedFiles

# Carichiamo le definizioni strutturali
include("source/NS2D_defs.jl")

# Prepariamo la memoria
global input_file = "U_channel.jl"
global info = Sinfo()
global first_run = true

# TRUCCO ANTI-WORLD AGE: Pre-carichiamo l'input nel Main.
# In questo modo variabili come 'file' e 'material' esistono già prima
# che NS2D.jl venga compilato, prevenendo il crash di World Age.
include("input/" * input_file)

# Ora possiamo avviare in sicurezza il solutore originale
include("NS2D.jl")
"""

# Scriviamo fisicamente il file launcher.jl
write("launcher.jl", launcher_code)
println("File launcher.jl preparato con successo.\n")


# =====================================================================
# CICLO PRINCIPALE DELLE SIMULAZIONI
# =====================================================================
for (i, params) in enumerate(Simulations)
    println("\n=======================================================")
    println(" SIMULAZIONE $i: Re = $(params.Re) | max_order = $(params.max_order)")
    println("=======================================================\n")

    # 1. LETTURA E MODIFICA DINAMICA DEL FILE
    origin = read(path_input, String)
    new = replace(origin, r"Re\s*=\s*[0-9\.]+" => "Re = $(params.Re)")
    new = replace(new, r"info\.max_order\s*=\s*[0-9]+" => "info.max_order = $(params.max_order)")
    new = replace(new, r"info\.Lmm\s*=\s*\[.*?\]" => "info.Lmm = $(params.Lmm)")
    new = replace(new, r"info\.Lmmcj\s*=\s*\[.*?\]" => "info.Lmmcj = $(params.Lmmcj)")
    write(path_input, new)
    println("File U_channel.jl aggiornato correttamente.")

    # 2. ESECUZIONE ISOLATA TRAMITE IL LAUNCHER
    try
        println("Avvio del calcolo...\n")
        
        # Lanciamo il nostro launcher preparatore
        run(`julia --project=. launcher.jl`) 
        
        println("\n ---> Simulazione $i Completata con successo!")

        # 3. SALVATAGGIO
        Re_clean = Int(params.Re)
        name = "results/Re$(Re_clean)_Ordine$(params.max_order)"
        
        println("Salvataggio in corso nella cartella: ", name)
        cp("output", name, force=true)

    catch e
        println("\n !!! ERRORE NELLA SIMULAZIONE $i !!!")
        println("La simulazione è fallita o è stata interrotta. Passo ai dati successivi.")
    end
end

# Pulizia finale
rm("launcher.jl", force=true)
println("\n*** TUTTE LE SIMULAZIONI IN CODA SONO TERMINATE ***")