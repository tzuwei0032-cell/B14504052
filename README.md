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
