import tkinter as tk
from tkinter import filedialog, messagebox
from PIL import Image, ImageTk
from sofistik_tools import SofiFileHandler, Ties, Plot, CDBData

def select_dat_file():
    """Ouvre une boîte de dialogue pour sélectionner un fichier .dat."""
    file_path = filedialog.askopenfilename(
        title="Sélectionner le fichier .dat",
        filetypes=[("Fichiers DAT", "*.dat"), ("Tous les fichiers", "*.*")]
    )
    dat_file_var.set(file_path)

def select_cdb_file():
    """Ouvre une boîte de dialogue pour sélectionner un fichier .cdb."""
    file_path = filedialog.askopenfilename(
        title="Sélectionner le fichier .cdb",
        filetypes=[("Fichiers CDB", "*.cdb"), ("Tous les fichiers", "*.*")]
    )
    cdb_file_var.set(file_path)

def on_next_click():
    """Fonction appelée lorsque l'utilisateur clique sur 'Suivant'."""
    dat_file = dat_file_var.get()
    cdb_file = cdb_file_var.get()
    if not dat_file or not cdb_file:
        messagebox.showwarning("Fichiers manquants", "Veuillez sélectionner les fichiers .dat et .cdb.")
        return
    
    try:
        # Créer une instance de SofiFileHandler pour manipuler le fichier .dat
        file_handler = SofiFileHandler()
        file_handler.add_dat(dat_file)

        # Ajouter un bloc de calcul linéaire
        n = 10  # Exemple de valeur pour LC
        file_handler.add_linear_calculation_block(n)
        file_handler.calculate_with_sps()  # Exécuter le calcul avec SOFiSTiK en mode batch

        # Extraire les données du CDB
        cdb_data = CDBData()
        cdb_data.open_cdb(cdb_file)
        node_coords = cdb_data.extract_node_coords()
        elements, element_numbers = cdb_data.extract_quad_elements()
        bounds = cdb_data.extract_structure_bounds()
        membran_forces = cdb_data.extract_quad_membran_forces()
        cdb_data.close_cdb()

        # Créer les tirants avec la classe Ties
        ties = Ties()
        threshold = 10  # Exemple de seuil de force
        generated_ties = ties.generate_reinforcement_ties_from_concrete(
            node_coords, elements, element_numbers, bounds, membran_forces, threshold, dir='z', n=10, m=10
        )

        # Afficher la structure avec les tirants conservés
        plot = Plot()
        fig = plot.plot_structure_with_ties(node_coords, elements, generated_ties, view_axis='xz')
        fig.show()  # Afficher la figure dans une nouvelle fenêtre

    except Exception as e:
        messagebox.showerror("Erreur d'exécution", f"Une erreur est survenue : {str(e)}")

# Configuration de la fenêtre principale
root = tk.Tk()
root.title("Interface de l'Application")
root.geometry("800x600")

# Charger l'image de fond
background_image_path = "images/pini_group_logo.jpg"  # Chemin de l'image
bg_image = Image.open(background_image_path)
bg_image = bg_image.resize((800, 600), Image.LANCZOS)
bg_photo = ImageTk.PhotoImage(bg_image)

# Créer un Canvas pour afficher l'image de fond
canvas = tk.Canvas(root, width=800, height=600)
canvas.pack(fill="both", expand=True)
canvas.create_image(0, 0, image=bg_photo, anchor="nw")

# Instructions pour l'utilisateur, centrées et espacées
instructions = tk.Label(canvas, text=(
    "1. Assurez-vous que votre structure est correctement configurée dans SOFiSTiK.\n"
    "2. Exportez votre fichier en .dat avec les matériaux, la structure SOFIMSHC et les charges SOFILOAD.\n"
    "3. Cliquez sur 'Suivant' pour sélectionner vos fichiers .dat et .cdb."
), justify="center", bg="white", font=("Arial", 12))
canvas.create_window(400, 80, window=instructions)

# Variables pour les chemins de fichiers
dat_file_var = tk.StringVar()
cdb_file_var = tk.StringVar()

# Bouton pour sélectionner le fichier .dat
dat_button = tk.Button(root, text="Sélectionner le fichier .dat", command=select_dat_file)
canvas.create_window(400, 220, window=dat_button)

# Label pour afficher le chemin du fichier .dat sélectionné
dat_label = tk.Label(root, textvariable=dat_file_var, fg="blue", bg="white")
canvas.create_window(400, 260, window=dat_label)

# Bouton pour sélectionner le fichier .cdb
cdb_button = tk.Button(root, text="Sélectionner le fichier .cdb", command=select_cdb_file)
canvas.create_window(400, 310, window=cdb_button)

# Label pour afficher le chemin du fichier .cdb sélectionné
cdb_label = tk.Label(root, textvariable=cdb_file_var, fg="blue", bg="white")
canvas.create_window(400, 350, window=cdb_label)

# Bouton 'Suivant' pour passer à l'étape suivante, centré en bas
next_button = tk.Button(root, text="Suivant", command=on_next_click)
canvas.create_window(400, 420, window=next_button)

# Lancer l'application Tkinter
root.mainloop()
