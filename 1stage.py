import numpy as np
import matplotlib.pyplot as plt
import random


# Физические константы
LIPOSOME_DIAMETER_NM = 200  # Диаметр липосомы в нанометрах
VOXEL_SIZE_NM = 10          # Размер одной ячейки (вокселя) в нанометрах
AVOGADRO_NUMBER = 6.022e23  # Число Авогадро

# Концентрации из системы PURE 
CONCENTRATIONS = {
    'RNA_POLYMERASE': 0.1,  # µM
    'RIBOSOME': 1.5,        # µM
    'ATP': 1500,            # µM (1.5 mM)
    'AMINO_ACID': 3000      # µM (3.0 mM)
}

grid_dim = int(LIPOSOME_DIAMETER_NM / VOXEL_SIZE_NM)
GRID_SIZE = (grid_dim, grid_dim, grid_dim)
volume_liters = (4/3) * np.pi * ((LIPOSOME_DIAMETER_NM / 2 * 1e-9)**3) * 1000

INITIAL_MOLECULES = {}
for name, conc_uM in CONCENTRATIONS.items():
    conc_mol_per_liter = conc_uM * 1e-6
    num_molecules = int(conc_mol_per_liter * volume_liters * AVOGADRO_NUMBER)
    INITIAL_MOLECULES[name] = num_molecules

print("--- Расчетные Начальные Условия ---")
print(f"Размер решетки: {GRID_SIZE}")
print(f"Объем липосомы (л): {volume_liters:.2e}")
print(f"Начальное количество молекул: {INITIAL_MOLECULES}")
print("---------------------------------")


EMPTY = 0
DNA_PROMOTER = 1
DNA_CODING = 2
RNA_POLYMERASE = 3
BOUND_RNAP = 4
RIBOSOME = 5
BOUND_RIBOSOME = 6
MRNA = 7
GFP_PROTEIN = 8
ATP = 9
AMINO_ACID = 10


class QuasiCellAutomaton:
    """Класс для моделирования синтеза белка в квази-клетке."""
    def __init__(self, grid_size=GRID_SIZE, params=None):
        self.grid_size = grid_size
        self.grid = np.zeros(self.grid_size, dtype=int)
        
        if params is None:
            self.params = {
                'p_diff': 0.2,
                'p_tscr_init': 0.1, 'p_tscr_elong': 0.9,
                'p_transl_init': 0.08, 'p_transl_elong': 0.8,
                'p_term': 0.01, 'p_atp_gen': 0.05,
                'p_mrna_degrad': 0.002, 'p_prot_degrad': 0.0001
            }
        else:
            self.params = params
            
        self.history = {'gfp': [], 'mrna': [], 'atp': []}
        self.initialize_grid()

    def initialize_grid(self):
        """Инициализация решетки с начальными молекулами."""
        x, y, z = [s // 2 for s in self.grid_size]
        self.grid[x, y, z] = DNA_PROMOTER
        for i in range(1, 6):
            if x + i < self.grid_size[0]:
                self.grid[x + i, y, z] = DNA_CODING

        empty_cells = list(np.argwhere(self.grid == EMPTY))
        random.shuffle(empty_cells) 

        for molecule_name, count in INITIAL_MOLECULES.items():
            molecule_id = globals()[molecule_name.upper()]
            placeable_count = min(count, len(empty_cells))
            for i in range(placeable_count):
                pos = tuple(empty_cells.pop())
                self.grid[pos] = molecule_id

    def get_neighbors(self, x, y, z):
        neighbors = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    if dx == 0 and dy == 0 and dz == 0:
                        continue
                    nx, ny, nz = x + dx, y + dy, z + dz
                    if 0 <= nx < self.grid_size[0] and 0 <= ny < self.grid_size[1] and 0 <= nz < self.grid_size[2]:
                        neighbors.append((nx, ny, nz))
        return neighbors

    def run_simulation(self, timesteps):
        """Запуск симуляции."""
        for t in range(timesteps):
            grid_to_read = np.copy(self.grid)
            
            cell_indices = list(np.ndindex(self.grid_size))
            random.shuffle(cell_indices)
            
            for x, y, z in cell_indices:

                if grid_to_read[x,y,z] == BOUND_RIBOSOME:

                    pass

                if grid_to_read[x,y,z] == BOUND_RNAP:

                    pass
                

            self.history['gfp'].append(np.sum(self.grid == GFP_PROTEIN))
            self.history['mrna'].append(np.sum(self.grid == MRNA))
            self.history['atp'].append(np.sum(self.grid == ATP))
            if t % 50 == 0:
                print(f"Шаг {t}: GFP={self.history['gfp'][-1]}, mRNA={self.history['mrna'][-1]}, ATP={self.history['atp'][-1]}")

    def plot_results(self):
        plt.figure(figsize=(12, 7))
        plt.plot(self.history['gfp'], label="GFP")
        plt.plot(self.history['mrna'], label="mRNA")
        plt.plot(self.history['atp'], label="ATP")
        plt.title("Динамика синтеза белка")
        plt.xlabel("Время")
        plt.ylabel("Количество")
        plt.legend()
        plt.grid(True)
        plt.show()

if __name__ == '__main__':
    sim = QuasiCellAutomaton()
    
    sim.run_simulation(2000)
    
    sim.plot_results()
