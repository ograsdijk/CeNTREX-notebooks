import papermill as pm
import centrex_TlF as centrex

fname_notebook_template = "multipass amherst template.ipynb"

transitions = [
    centrex.transitions.LaserTransition("R1", F1=3/2, F=1),
    centrex.transitions.LaserTransition("R1", F1=3/2, F=2),
    centrex.transitions.LaserTransition("R1", F1=5/2, F=2),
    centrex.transitions.LaserTransition("R1", F1=5/2, F=3),
    centrex.transitions.LaserTransition("Q1", F1=1/2, F=0),
    centrex.transitions.LaserTransition("Q1", F1=1/2, F=1),
    centrex.transitions.LaserTransition("Q1", F1=3/2, F=1),
]

Js_in_systems = [
    [1,3],
    [1,3],
    [1,3],
    [1,3],
    [1,3],
    [1,3],
    [1,3],
]

for transition, Js_in_system in zip(transitions, Js_in_systems):
    print(f"Executing notebook for {transition}")
    fname_notebook = f"{transition}"
    fname_notebook = fname_notebook.replace(",","").replace("/","_")
    fname_notebook = fname_notebook.replace("Transition(", "")
    fname_notebook = fname_notebook[:-1]
    fname_notebook = f"{fname_notebook}.ipynb"
    
    pm.execute_notebook(
        fname_notebook_template,
        fname_notebook,
        parameters = dict(
            transition_type = transition.transition,
            excited_J = transition.excited_selector.J,
            excited_F = transition.excited_selector.F,
            excited_F1 = transition.excited_selector.F1,
            system_Js = Js_in_system,
        )
    )
    