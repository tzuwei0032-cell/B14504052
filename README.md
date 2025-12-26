程式功能與技術原理
	本程式功能在於計算酸鹼滴定過程及到達低定點後，容液中各離子濃度的實際值(不經過取近似)，以溶液店中性理論發展成一套公式，並生成pH值滴定曲線圖。

	
使用方法
1.	啟動程式
2.	選擇酸鹼種類的輸入方式(手動/資料庫)
3.	依照指示輸入相關資訊
4.	輸入被滴定液初始濃度
5.	輸入被滴定液初始體積
6.	輸入滴定液初始濃度
7.	欲討論童狀況下不同滴定液體積之情形則接續照指示輸入滴定液體積即可，不須重啟城市；若欲察看pH值曲線圖形則輸入q
8.	輸入q後即生成pH直曲線圖


程式架構



開發過程
    使用者輸入被滴定液的相關數據c0(初始濃度)、v0(初始體積)、ks(溶液的各級平衡常數；判斷酸鹼元數－n=輸入的K值數量→
    處理K值，配合結果以pH值呈現，並方便統一計算，先判斷被滴定物的酸鹼，將所有K值轉為Ka值(※這邊出錯了－Kb的轉換應用reversed()因為Kb1對應的是解離常數的最後一階) 
	
生成離子之化學式→

進入離子濃度計算－先建構公式→

建構公式→
	首先，我們回顧過去我們曾經學過二、三原酸的公式解解法然而其計算過程使用了近似，因此不適用於本term projec；或許，我們在不用近似的情況下解出一元二次、三次方程式的解亦可達成目標，然而此法並不具備通用性，因此我們希望尋找一個具備通用性的解法滿足1~N元酸的酸解滴定。
	電荷平衡 : 溶液達電荷平衡
		[H^+ ]+[〖Na〗^+ ]=[〖OH〗^- ]+∑▒〖i*[A^(i-) ] 〗
			質量守恆：所有解離離子的濃度總和等於初始濃度
			分率函數：定義α_i為第i級解離物種佔總濃度的比例。這是一個依賴[H^+ ]的複雜多項式。
			以三元酸為例
			各離子濃度
		K_a1=[H^+ ][H_2 A^- ]/[H_3 A]   → [H_2 A^- ]=K_a1/[H^+ ]  [H_3 A]
		K_a2=[H^+ ][HA^(2-) ]/[H_2 A^- ]   → [HA^(2-) ]=(K_a1 K_a2)/[H^+ ]^2  [H_3 A]
		K_a3=[H^+ ][A^(3-) ]/[HA^(2-) ]   → [A^(3-) ]=(K_a1 K_a2 K_a3)/[H^+ ]^3  [H_3 A]
		C_t=[H_3 A]+[H_2 A^- ]+[HA^(2-) ]+[A^(3-) ]
			=[H_3 A](1+K_a1/[H^+ ] +(K_a1 K_a2)/[H^+ ]^2 +(K_a1 K_a2 K_a3)/[H^+ ]^3 )
			令D=[H^+ ]^3+[H^+ ]^2 K_a1+[H^+ ] K_a1 K_a2+K_a1 K_a2 K_a3
		各離子的表達式
			未解離α_0=[H^+ ]/D
			一級解離α_1=([H^+ ]^2 K_a1)/D
			二級解離α_2=([H^+ ] K_a1 K_a2)/D
			三級解離α_2=(K_a1 K_a2 K_a3)/D

		利用電荷平衡方程式尋找[H^+ ]→

		先求得混合後濃度→

	電荷平衡→

	酸被強鹼滴定：[H+] + [Na+] = [OH-] + Σ(i * [A^i-])→

	鹼被強酸滴定：[H+] + Σ((n-i) * [BH+^i]) = [OH-] + [Cl-]→

	利用brentq求出[H+]→

	輸出→

	使用者輸入→

	重複滴定，不同滴定體積(相同情況下)→

參考資料來源
	1. 使用LLM詢問離子濃度部曲近似的計算方法
	2. Kb轉Ka時的除錯
	3. brentq的用法

增強與改良
	最初的半本並未有資料庫功能與繪圖功能，先是在做出第一版後事跑時注意到，要使用本程式總是要先上網搜尋K值，然而我們在做實驗時很多時候手邊並未有準確的K值數據，於是建構一個裝有大量酸鹼之K值的資料夠將能使此程式更方便
	又在將入資料庫後發現似乎可以再加入pH值曲線圖的繪圖功能來完善此程式，使得此程式更加視覺化
	再加入圖表後又發現，應增加滴定點後過量滴定液的樣本以完善圖表
	








Three introduction of the code
This code include three core parts:

1. Chemical Database
No more flipping through textbooks to find dissociation constants(K). The program features a built-in library of common acids and bases. Simply enter the name (e.g., Acetic acid or Phosphoric acid), and the system automatically configures the chemical properties for you.

2. High-Precision pH Solver
Instead of using simplified approximation formulas, this simulator employs the Charge Balance Method. By solving complex equilibrium equations numerically, it provides highly accurate pH values regardless of how dilute, strong, or weak the solution is.
,  Dynamic Titration Curve Visualization

3.The program automatically simulates 300 data points to generate a smooth pH curve. It strategically marks the "Equivalence Points" with red dots, allowing you to visualize exactly when the chemical neutralization is complete.


How to Operate (3 Simple Steps)

Step 1: Define Your "Analyte"
When the program starts, choose your input mode:
	Easy Mode: Select 1 and enter the chemical name (e.g., Acetic acid).
	Custom Mode: Select 2 to manually input the $K$ values (dissociation constants).
    
Step 2: Set Experimental Parameters
Think of this as preparing your lab equipment. You will need to provide:
	Concentration (C_0): How concentrated is your sample (in Molarity)?
	Initial Volume (V_0): How many mL are in your beaker?
	Titrant Concentration (C_t): The concentration of the strong acid or base in the burette.
    
Step 3: Analyze Data & View Plot
	Query Data: Enter any volume (mL) to see the instantaneous pH and a detailed breakdown of all chemical species in the solution.
	Generate Plot: Type q to finish. A high-quality titration graph will pop up, showing you exactly where the pH "jumps" occur.
