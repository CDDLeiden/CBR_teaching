import time
import uuid
import json

import pandas as pd
import pyperclip


def show_residues_around(view, component_index=0, selection="ligand", radius=5.0):
    js = (
        f"""
        // Get first (and only) loaded component: our protein-ligand system
        var system = this.stage.compList[{component_index}]; 
        // Store current selection, we will need it later
        var prevSele = system.selection.string;
        // Set selection to our desired ligand
        system.setSelection("{selection}");
        // Select all atoms within 5Å from the ligand
        var around = system.structure.getAtomSetWithinSelection(system.selection, {radius});
        """
        # Not f-string below
        """
        // Extend selection so it includes full residues
        var around_complete = system.structure.getAtomSetWithinGroup(around);
        // Add representation for those atoms
        system.addRepresentation("licorice", {sele: around_complete.toSeleString()});
        // Restore selection to original one; otherwise changes won't be reflected
        system.setSelection(prevSele)
        """
    )
    view._execute_js_code(js)



def get_residues_around(view, component_index=0, selection="ligand", radius=5.0, timeout: int = 20):
    """
    Executes JavaScript to calculate surrounding residues and copy them to the
    system clipboard. This function DOES NOT wait for a return value.

    :param timeout: Max seconds to wait for the JavaScript to finish.
    """
    # Generate a unique token to signal that we are waiting for a result.
    wait_token = f"nglview-waiting-{uuid.uuid4().hex}"
    try:
        pyperclip.copy(wait_token)
    except pyperclip.PyperclipException as e:
        print("PYPERCLIP ERROR: Could not access the system clipboard.")
        print("Please ensure you have a clipboard utility installed (e.g., xclip on Linux).")
        print(f"Error details: {e}")
        return []

    js = (
        f"""
        function runWhenReady() {{
            // Get first (and only) loaded component: our protein-ligand system
            var system = this.stage.compList[{component_index}];
            
            if (system && system.structure && system.structure.atomCount > 0) {{
                // Store current selection, we will need it later
                var prevSele = system.selection.string;
                // Set selection to our desired ligand
                system.setSelection("{selection}");
                // Select all atoms within 5Å from the ligand
                var around = system.structure.getAtomSetWithinSelection(system.selection, {radius});
                """
                # Not f-string below
                """
                // Extend selection so it includes full residues
                var around_complete = system.structure.getAtomSetWithinGroup(around);
                var around_complete_atoms = [];
                around_complete.forEach(item => around_complete_atoms.push(item));
                // Restore selection to original one; otherwise changes won't be reflected
                system.setSelection(prevSele)
                // Store residue data
                var residueData = [];
                var uniqueResidues = new Set();
                system.structure.eachAtom(atom => {
                                                var residue = atom.residue;
                                                if (atom && residue && around_complete_atoms.includes(atom.serial)) {
                                                    var residueKey = residue.chainname + ':' + residue.resno;
                                                    if (!uniqueResidues.has(residueKey)) {
                                                        uniqueResidues.add(residueKey);
                                                        residueData.push({
                                                            'resi': residue.resno,
                                                            'resn': residue.resname,
                                                            'chain': residue.chainname
                                                        });
                                                    }
                                                }
                                            });
                var jsonData = JSON.stringify(residueData);
                navigator.clipboard.writeText(jsonData);
            } else {
                setTimeout(runWhenReady.bind(this), 200);
            }
        }
        runWhenReady.bind(this)();
        """
    )

    # Step 2: Execute the JavaScript and start polling
    view._execute_js_code(js)
    start_time = time.time()
    while True:
        # Check if the clipboard has been updated
        clipboard_content = pyperclip.paste()
        if clipboard_content != wait_token:
            break
        # Check for timeout
        if time.time() - start_time > timeout:
            print(f"ERROR: Timed out after {timeout} seconds.")
            return []
        # Wait a short moment before polling again to avoid excessive CPU usage
        time.sleep(0.5)
    # Step 3: Parse and return the result
    try:
        residue_list = json.loads(clipboard_content)
        return pd.DataFrame(residue_list)
    except (json.JSONDecodeError, TypeError):
        print(f"ERROR: JavaScript encountered an error.")
        return []
