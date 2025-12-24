import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

CHEMICAL_DATABASE = {
    "Hydrochloric acid": {"type": "A", "ks": [1.3e6]},
    "Nitric acid": {"type": "A", "ks": [2.4e1]},
    "Perchloric acid": {"type": "A", "ks": [1e9]},
    "Sulfuric acid": {"type": "A", "ks": [1e3, 1.2e-2]},

    "Phosphoric acid": {"type": "A", "ks": [7.11e-3, 6.32e-8, 4.22e-13]},
    "Carbonic acid": {"type": "A", "ks": [4.3e-7, 5.6e-11]},
    "Boric acid": {"type": "A", "ks": [5.8e-10]},
    "Hydrofluoric acid": {"type": "A", "ks": [6.6e-4]},
    "Hydrocyanic acid": {"type": "A", "ks": [6.2e-10]},
    "Nitrous acid": {"type": "A", "ks": [7.2e-4]},
    "Sulfurous acid": {"type": "A", "ks": [1.3e-2, 6.2e-8]},
    "Hypochlorous acid": {"type": "A", "ks": [2.9e-8]},

    "Acetic acid": {"type": "A", "ks": [1.75e-5]},
    "Formic acid": {"type": "A", "ks": [1.8e-4]},
    "Oxalic acid": {"type": "A", "ks": [5.6e-2, 5.42e-5]},
    "Citric acid": {"type": "A", "ks": [7.4e-4, 1.7e-5, 4.0e-7]},
    "Benzoic acid": {"type": "A", "ks": [6.3e-5]},
    "Lactic acid": {"type": "A", "ks": [1.38e-4]},
    "Ascorbic acid": {"type": "A", "ks": [8.0e-5, 1.6e-12]},

    "Ammonia": {"type": "B", "ks": [1.8e-5]},
    "Ethylenediamine": {"type": "B", "ks": [8.5e-5, 7.1e-8]},
    "Pyridine": {"type": "B", "ks": [1.8e-9]},
    "Aniline": {"type": "B", "ks": [4.3e-10]},
    "Methylamine": {"type": "B", "ks": [4.4e-4]},
    "Trimethylamine": {"type": "B", "ks": [6.4e-5]},
    
    "Sodium hydroxide": {"type": "B", "ks": [1e9]},
    "Potassium hydroxide": {"type": "B", "ks": [1e9]}
}

class UniversalTitrator:
    def __init__(self, c0, v0, ks, is_acid=True):
        self.c0 = c0
        self.v0 = v0
        self.kw = 1e-14
        self.is_acid = is_acid
        n = len(ks)

        if is_acid:
            self.kas = ks
        else:
            self.kas = [(self.kw / kb) for kb in reversed(ks)]
        
        if is_acid:
            self.labels = [f"H{n-i}A{'' if i==0 else '-'*i}" for i in range(n+1)]
        else:
            self.labels = [f"BH{n-i}{'' if i==n else '+'*(n-i)}" for i in range(n+1)]

    def get_alphas(self, h_conc):
        n = len(self.kas)
        coeffs = [1.0]
        temp = 1.0
        for ka in self.kas:
            temp *= ka
            coeffs.append(temp)
        
        terms = [coeffs[i] * (h_conc**(n - i)) for i in range(n + 1)]
        denominator = sum(terms)
        return [t / denominator for t in terms]

    def solve_ph(self, v_titrant, c_titrant):
        v_total = self.v0 + v_titrant
        curr_c_sub = self.c0 * self.v0 / v_total
        curr_c_tit = c_titrant * v_titrant / v_total

        def charge_balance_error(h_conc):
            oh = self.kw / h_conc
            alphas = self.get_alphas(h_conc)
            if self.is_acid:
                return (h_conc + curr_c_tit) - (oh + sum(i * alphas[i] for i in range(1, len(alphas))) * curr_c_sub)
            else:
                return (h_conc + sum((len(alphas)-1-i) * alphas[i] for i in range(len(alphas)-1)) * curr_c_sub) - (oh + curr_c_tit)

        return brentq(charge_balance_error, 1e-16, 15.0)

def main():
    print("="*50)
    print("      多元酸鹼精確滴定求解器 (修正版)")
    print("="*50)
    
    try:
        print("\n[選項] 1. 使用資料庫 2. 手動輸入")
        choice = input("請選擇 (1/2): ")
        
        if choice == '1':
            name = input("請輸入物質名稱 (如 Acetic acid): ")
            if name in CHEMICAL_DATABASE:
                entry = CHEMICAL_DATABASE[name]
                is_acid = (entry["type"] == 'A')
                ks = entry["ks"]
            else:
                print("找不到該物質！")
                return
        else:
            sub_type = input("滴定對象是酸 (A) 還是鹼 (B)? ").upper()
            is_acid = (sub_type == 'A')
            k_input = input("請輸入各級 K 值 (以逗號隔開): ")
            ks = [float(k.strip()) for k in k_input.split(',')]

        c0 = float(input("滴定對象初始濃度 (M): "))
        v0 = float(input("滴定對象初始體積 (mL): "))
        ct = float(input("滴定液 (強酸/強鹼) 濃度 (M): "))
        
        titrator = UniversalTitrator(c0, v0, ks, is_acid)
        n_protic = len(ks)
        v_eq_final = (c0 * v0 * n_protic / ct) * 1.5
        v_plot = np.linspace(0, v_eq_final, 300) 
        ph_plot = []
        
        print("\n系統正在計算繪圖點位...")
        for vp in v_plot:
            h_sol = titrator.solve_ph(vp, ct)
            ph_plot.append(-np.log10(h_sol))
        plt.figure(figsize=(10, 6))
        plt.plot(v_plot, ph_plot, label='Titration Curve', color='blue', lw=2)
        for i in range(1, n_protic + 1):
            v_eq = (c0 * v0 * i) / ct
            h_eq = titrator.solve_ph(v_eq, ct)
            ph_eq = -np.log10(h_eq)
            plt.plot(v_eq, ph_eq, 'ro') 
            plt.annotate(f'Eq {i}\n({v_eq:.1f}mL, pH{ph_eq:.1f})', 
                         xy=(v_eq, ph_eq), xytext=(v_eq + v_eq_final*0.05, ph_eq),
                         arrowprops=dict(arrowstyle='->', color='red'))

        plt.title(f"Titration Curve of {'Acid' if is_acid else 'Base'}")
        plt.xlabel("Volume of Titrant Added (mL)")
        plt.ylabel("pH")
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.legend()
        print(" 完成繪圖 ")
        v_equivalence = (c0 * v0 * len(ks)) / ct
        print(f"\n[資訊] 理論完全當量點體積: {v_equivalence:.2f} mL")
        
        while True:
            target_v = input("\n輸入滴定液體積 (mL)，或按 'q' 顯示圖表並結束: ")
            if target_v.lower() == 'q': break
            
            v = float(target_v)
            h_final = titrator.solve_ph(v, ct)
            alphas = titrator.get_alphas(h_final)
            curr_total_c = (c0 * v0) / (v0 + v)
            
            print(f"\n--- 結果報告 (加入 {v} mL) ---")
            print(f"{'物種':<18} | {'濃度 (M)':<15}")
            print("-" * 35)
            print(f"{'[H+]':<18} | {h_final:.4e}")
            print(f"{'[OH-]':<18} | {(1e-14/h_final):.4e}")
            
            for i, label in enumerate(titrator.labels):
                conc = curr_total_c * alphas[i]
                print(f"{'['+label+']':<18} | {conc:.4e}")
            
            print("-" * 35)
            print(f"最終 pH 值: {-np.log10(h_final):.4f}")
        plt.show()

    except Exception as e:
        print(f"執行出錯: {e}")


if __name__ == "__main__":
    main()
