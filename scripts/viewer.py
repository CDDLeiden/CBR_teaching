import io
import time
import uuid
import json
import threading
import itertools

import pandas as pd
#import pyperclip
import nglview as nv
from IPython.display import display, HTML
import Bio.PDB


def show_residues_around(view, component_index=0, selection="ligand", radius=5.0):
    js = (
        f"""
        // Get first (and only) loaded component: our protein-ligand system
        var system = this.stage.compList[{component_index}]; 
        // Store current selection, we will need it later
        var prevSele = system.selection.string;
        // Set selection to our desired ligand
        system.setSelection("{selection}");
        // Select all atoms within 5A from the ligand
        var around = system.structure.getAtomSetWithinSelection(system.selection, {radius});
        """
        """
        // Extend selection so it includes full residues
        var around_complete = system.structure.getAtomSetWithinGroup(around);
        // Add representation for those atoms
        system.addRepresentation("licorice", {sele: around_complete.toSeleString()});
        // Restore selection to the original one; otherwise, changes won't be reflected
        system.setSelection(prevSele)
        """
    )
    view._execute_js_code(js)

def show_residues(view, component_index=0, selection="protein"):
    js = (
        f"""
        // Get first (and only) loaded component: our protein-ligand system
        var system = this.stage.compList[{component_index}]; 
        // Store current selection, we will need it later
        var prevSele = system.selection.string;
        // Set selection to our desired ligand
        system.setSelection("{selection}");
        """
        """
        // Add representation for those atoms
        system.addRepresentation("licorice", {sele: system.selection.toSeleString()});
        // Restore selection to the original one; otherwise, changes won't be reflected
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


def show_residues_around(view, component_index=0, selection="ligand", radius=5.0):
    js = (
        f"""
        // Get first (and only) loaded component: our protein-ligand system
        var system = this.stage.compList[{component_index}]; 
        // Store current selection, we will need itlaterr
        var prevSele = system.selection.string;
        // Set selection to our desired ligand
        system.setSelection("{selection}");
        // Select all atoms within 5A from the ligand
        var around = system.structure.getAtomSetWithinSelection(system.selection, {radius});
        """
        """
        // Extend selection so it includes full residues
        var around_complete = system.structure.getAtomSetWithinGroup(around);
        // Add representation for those atoms
        system.addRepresentation("licorice", {sele: around_complete.toSeleString()});
        // Restore selection to the original one; otherwise, changes won't be reflected
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


class NGLViewResidueHelper:
    """
    A helper class that provides a robust method to get data back from NGLView.

    It works by establishing a persistent message handler, which solves race
    conditions present in some Jupyter environments where messages from asynchronous
    JavaScript calls can be dropped.

    Workflow:
    1. Create an instance: `helper = NGLViewResidueHelper(view)`
    2. Call the method: `residues = helper.get_residues_around(selection="UNL")`
    """
    def __init__(self, view):
        self.view = view
        self._results = {}
        # Establish a single, persistent message handler for this view object.
        self.view.on_msg(self._on_msg)

    def _on_msg(self, widget, msg, buffers):
        """This function runs every time ANY message comes from the JS view."""
        # Check if it's the specific type of message we're looking for
        if msg.get("type") == "ngl_residue_result":
            request_id = msg.get("request_id")
            if request_id in self._results:
                # Place the data in the correct slot and signal that it has arrived.
                self._results[request_id]["data"] = msg.get("data", [])
                self._results[request_id]["event"].set()

    def get_residues_around(self, component_index=0, selection="ligand", radius=5.0, timeout=20):
        """
        Calculates and returns a list of residues surrounding a selection in a
        single, blocking call using a robust internal communication channel.

        :param selection: The residue name to select (e.g., "UNL", "ALA", "ligand").
        :param radius: The radius in Ångströms to select residues within.
        :param timeout: Max seconds to wait for the JavaScript to return a result.
        :return: A list of dictionaries representing surrounding residues.
        """
        # Generate a unique ID for this specific request.
        request_id = f"req_{uuid.uuid4().hex}"
        event = threading.Event()
        self._results[request_id] = {"event": event, "data": []}
        # Format the selection string robustly.
        if not selection.startswith(":") and selection.upper() != 'LIGAND':
            formatted_selection = ":" + selection
        else:
            formatted_selection = selection
        # This is YOUR "bulletproof" JavaScript logic, with one addition:
        # it now sends the unique request_id back with the data.
        js = f"""
        function runWhenReady() {{
            var system = this.stage.compList[{component_index}];
            if (system && system.structure && system.structure.atomCount > 0) {{
                console.log("SUCCESS: Structure loaded. Running selection logic...");
                try {{
                    system.setSelection('{formatted_selection}');
                    var around = system.structure.getAtomSetWithinSelection(system.selection, {radius});
                    var around_complete = system.structure.getAtomSetWithinGroup(around);
                    
                    var surrounding_atom_indices = [];
                    around_complete.forEach(atom_serial => surrounding_atom_indices.push(atom_serial));
                    
                    var residueData = [];
                    var uniqueResidues = new Set();
                    var targetResname = "{formatted_selection.replace(":", "")}";

                    system.structure.eachAtom(atom => {{
                        if (atom && atom.residue && surrounding_atom_indices.includes(atom.index)) {{
                            var residue = atom.residue;
                            if (residue.resname !== targetResname) {{
                                var residueKey = residue.chainname + ':' + residue.resno;
                                if (!uniqueResidues.has(residueKey)) {{
                                    uniqueResidues.add(residueKey);
                                    residueData.push({{
                                        'resi': residue.resno, 'resn': residue.resname, 'chain': residue.chainname
                                    }});
                                }}
                            }}
                        }}
                    }});
                    
                    // --- ROBUST COMMUNICATION ---
                    // Send the data back with the unique ID for this request.
                    console.log("SUCCESS: Calculation finished. Sending " + residueData.length + " residues back to Python.");
                    this.model.send({{
                        type: 'ngl_residue_result',
                        request_id: '{request_id}',
                        data: residueData
                    }});
                }} catch (e) {{
                    console.error("A FATAL ERROR OCCURRED IN THE JS CODE:", e);
                    this.model.send({{ type: 'ngl_residue_result', request_id: '{request_id}', data: ["ERROR"] }});
                }}
            }} else {{
                console.log("Waiting for structure to load...");
                setTimeout(runWhenReady.bind(this), 200);
            }}
        }}
        runWhenReady.bind(this)();
        """
        print(f"DEBUG (Python): Executing JS for selection '{formatted_selection}'. Waiting for result...")
        self.view._execute_js_code(js)
        # Wait for the _on_msg handler to signal that this specific request is done.
        received = event.wait(timeout=timeout)
        if not received:
            print(f"ERROR: Timed out after {timeout} seconds. No result received from JavaScript.")
            del self._results[request_id]
            return []
        # Retrieve the data, clean up, and return.
        result_data = self._results[request_id]["data"]
        del self._results[request_id]
        if result_data == ["ERROR"]:
            print("ERROR: The JavaScript code encountered a fatal error. Check the browser console.")
            return []
        print("Successfully received data from JavaScript.")
        return result_data


def get_residues_around(pdb_file, selection="LIG", chain=None, radius=5.0):
    """Determine residues surrounding a selection using BioPython.

    :param pdb_file: PDB file.
    :param selection: The 3-letter residue name of the ligand/residue to search around.
    :param chain: (Optional) The specific chain ID of the ligand, if needed.
    :param radius: The radius in Ångströms for the search.
    :return: A list of dictionaries, where each dictionary represents a surrounding residue.
             Returns None if the selection is not found.
    """
    # --- Step 1: Load the structure into BioPython ---
    pdb_text = ""
    # Check if the source is a PDB ID (4 characters)
    # print(f"Reading from file: {pdb_file}...")
    try:
        with open(pdb_file, 'r') as f:
            pdb_text = f.read()
    except FileNotFoundError:
        raise RuntimeError(f"Error: File not found at {pdb_file}")
    
    if not pdb_text:
        raise RuntimeError("Error: Could not load PDB data.")
    if pdb_file.lower().endswith('pdb') or pdb_file.lower.endswith('ent'):
        pdb_parser = Bio.PDB.PDBParser(QUIET=True)
        structure = pdb_parser.get_structure("my_structure", io.StringIO(pdb_text))
        model = structure[0]
    elif pdb_file.lower().endswith('cif'):
        cif_parser = Bio.PDB.FastMMCIFParser(QUIET=True)
        structure = cif_parser.get_structure("my_structure", io.StringIO(pdb_text))
        model = structure[0]
    else:
        raise RuntimeError("Error: Could not load PDB data.")

    # --- Step 2: Find the ligand and protein atoms ---
    atoms = []
    ligand_atoms = []
    protein_atoms = []
    
    for res in model.get_residues():
        # Check for the ligand
        if res.get_resname() == selection:
            if chain is None or res.get_parent().id == chain:
                ligand_atoms.extend(res.get_atoms())
        # Check for standard protein residues
        elif Bio.PDB.is_aa(res):
            protein_atoms.extend(res.get_atoms())
    atoms = ligand_atoms + protein_atoms

    if not ligand_atoms:
        raise RuntimeError(f"Error: Selection '{selection}' (Chain: {chain or 'any'}) not found in the structure.")
    
    # print(f"Found {len(ligand_atoms)} atoms for selection '{selection}'.")
    # print(f"Found {len(protein_atoms)} atoms for the protein.")
    
    # --- Step 3: Perform the neighbor search ---
    ns = Bio.PDB.NeighborSearch(protein_atoms)
    # Find all protein atoms within `radius` of any ligand atom
    nearby_atoms = set()
    for atom in ligand_atoms:
        neighbors = ns.search(atom.coord, radius, level='A') # 'A' for atom level
        nearby_atoms.update(neighbors)

    nearby_residues = set(atom.get_parent() for atom in nearby_atoms)

    # --- Step 4: Get the unique residues from the nearby atoms ---
    surrounding_residues = []
    for residue in nearby_residues:
        surrounding_residues.append({"resid": residue.get_id()[1],
                                     "resname": residue.get_resname(),
                                     "chain": residue.get_parent().id,
                                     })
    # print(f"Found {len(surrounding_residues)} surrounding residues.")

    return pd.DataFrame(surrounding_residues)