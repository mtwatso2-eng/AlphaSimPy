import pandas as pd
from smolagents import CodeAgent, load_tool, tool
import os
import json
import base64
from typing import Any, Optional
from smolagents.tools import Tool
from markdownify import markdownify
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
import requests


class LoadBraidProgramTool(Tool):
    """Loads a BRAID breeding program abstraction from disk."""
    name = "load_braid_program"
    description = "Loads a BRAID breeding program abstraction file from a provided path and returns its contents as text. Use this first when a braid_program_path is provided."
    inputs = {
        'braid_program_path': {
            'type': 'string',
            'description': 'Path to a BRAID abstraction file (YAML or JSON).'
        }
    }
    output_type = "string"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_initialized = False

    def forward(self, braid_program_path: str) -> str:
        if not os.path.exists(braid_program_path):
            return f"Error: BRAID program file not found at {braid_program_path}"
        try:
            with open(braid_program_path, "r", encoding="utf-8") as f:
                return f.read()
        except Exception as e:
            return f"Error loading BRAID program: {type(e).__name__}: {e}"


class SaveNotebookTool(Tool):
    """Saves a Python notebook (.ipynb) to disk."""
    name = "save_notebook"
    description = "Saves AlphaSimPy code as a Jupyter notebook (.ipynb file). Use this to write the final output notebook. Takes a list of dicts with 'cell_type' and 'source' keys."
    inputs = {
        'cells': {
            'type': 'object',
            'description': 'List of notebook cells. Each cell: {"cell_type": "code"|"markdown", "source": "string or list of strings"}'
        },
        'output_path': {
            'type': 'string',
            'description': 'Path where to save the .ipynb file'
        }
    }
    output_type = "string"

    def forward(self, cells: list, output_path: str) -> str:
        import nbformat as nbf
        nb = nbf.v4.new_notebook()
        for c in cells:
            src = c.get("source", "")
            if isinstance(src, list):
                src = "".join(src)
            ct = c.get("cell_type", "code")
            if ct == "markdown":
                nb.cells.append(nbf.v4.new_markdown_cell(src))
            else:
                nb.cells.append(nbf.v4.new_code_cell(src))
        with open(output_path, "w") as f:
            nbf.write(nb, f)
        return f"Notebook saved to {output_path}"


class SaveAlphaSimPyNotebookTool(Tool):
    """Saves an AlphaSimPy notebook (.ipynb) to disk."""
    name = "save_alphasimpy_notebook"
    description = (
        "Saves an AlphaSimPy Jupyter notebook (.ipynb). "
        "Pass notebook cells as a list of dicts with `cell_type` and `source`. "
        "Creates parent directories for output_path if needed. Returns the saved file path."
    )
    inputs = {
        'cells': {
            'type': 'object',
            'description': 'List of notebook cells. Each cell: {"cell_type": "code"|"markdown", "source": "string or list of strings"}'
        },
        'output_path': {
            'type': 'string',
            'description': 'Path where to save the AlphaSimPy notebook (.ipynb)'
        }
    }
    output_type = "string"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_initialized = False

    def forward(self, cells: list, output_path: str) -> str:
        import nbformat as nbf
        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        nb = nbf.v4.new_notebook()
        for c in cells:
            src = c.get("source", "")
            if isinstance(src, list):
                src = "".join(src)
            ct = c.get("cell_type", "code")
            if ct == "markdown":
                nb.cells.append(nbf.v4.new_markdown_cell(src))
            else:
                nb.cells.append(nbf.v4.new_code_cell(src))
        with open(output_path, "w", encoding="utf-8") as f:
            nbf.write(nb, f)
        return f"AlphaSimPy notebook saved to {output_path}"


class ListFolderContentsTool(Tool):
    name = "list_folder_contents"
    description = "Lists the full paths of all files and subfolders within a given folder path."
    inputs = {
        'folder_path': {
            'type': 'string',
            'description': 'Path to the folder.'
        }
    }
    output_type = "object"

    def forward(self, folder_path: str) -> object:
        all_paths = []
        for root, dirs, files in os.walk(folder_path):
            for d in dirs:
                all_paths.append(os.path.join(root, d))
            for f in files:
                all_paths.append(os.path.join(root, f))
        return all_paths

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_initialized = False


class ReadTextFileTool(Tool):
    name = "read_text_file"
    description = (
        "Reads a UTF-8 text file and returns its contents. "
        "Use this instead of Python open() in agent-generated code."
    )
    inputs = {
        'file_path': {
            'type': 'string',
            'description': 'Path to the text file to read.'
        }
    }
    output_type = "string"

    def forward(self, file_path: str) -> str:
        if not os.path.exists(file_path):
            return f"Error: File not found at {file_path}"
        if not os.path.isfile(file_path):
            return f"Error: Path is not a file: {file_path}"
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                return f.read()
        except Exception as e:
            return f"Error reading file: {type(e).__name__}: {e}"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_initialized = False

class FinalAnswerTool(Tool):
    name = "final_answer"
    description = "Provides a final answer to the given problem."
    inputs = {'answer': {'type': 'string', 'description': 'The final answer to the problem'}}
    output_type = "string"
    def forward(self, answer: str) -> str:
        return answer
    def __init__(self, *args, **kwargs):
        self.is_initialized = False

class RenderMermaidDiagramTool(Tool):
    name = "render_mermaid_diagram"
    description = "Renders a mermaid diagram code into a PNG image file. Takes mermaid diagram code as input and saves it as a PNG image. Returns the path to the saved image file."
    inputs = {
        'mermaid_code': {
            'type': 'string',
            'description': 'The mermaid diagram code to render (without the ```mermaid wrapper, just the diagram code itself)'
        },
        'output_path': {
            'type': 'string',
            'description': 'Optional. The path where to save the rendered diagram. If not provided, will save as breeding_program_diagram.png in the current directory.',
            'nullable': True
        }
    }
    output_type = "string"
    
    def forward(self, mermaid_code: str, output_path: Optional[str] = None) -> str:
        """
        Renders mermaid diagram code to a PNG image using the mermaid.ink API.
        """
        if output_path is None:
            output_path = "breeding_program_diagram.png"
        
        # Clean the mermaid code - remove any markdown code block wrappers if present
        mermaid_code = mermaid_code.strip()
        if mermaid_code.startswith("```mermaid"):
            mermaid_code = mermaid_code[9:].strip()
        if mermaid_code.startswith("```"):
            mermaid_code = mermaid_code[3:].strip()
        if mermaid_code.endswith("```"):
            mermaid_code = mermaid_code[:-3].strip()
        
        # Encode the mermaid code to base64
        mermaid_bytes = mermaid_code.encode('utf-8')
        mermaid_base64 = base64.urlsafe_b64encode(mermaid_bytes).decode('utf-8')
        
        # Use mermaid.ink API to render the diagram
        api_url = f"https://mermaid.ink/img/{mermaid_base64}"
        
        try:
            # Download the rendered image
            response = requests.get(api_url, timeout=30)
            response.raise_for_status()
            
            # Save the image
            with open(output_path, 'wb') as f:
                f.write(response.content)
            
            return f"Mermaid diagram rendered and saved to: {output_path}"
        except Exception as e:
            # Fallback: try using playwright if available
            try:
                from playwright.sync_api import sync_playwright
                
                html_content = f"""
                <!DOCTYPE html>
                <html>
                <head>
                    <script src="https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid.min.js"></script>
                    <script>
                        mermaid.initialize({{ startOnLoad: true }});
                    </script>
                </head>
                <body>
                    <div class="mermaid">
                    {mermaid_code}
                    </div>
                </body>
                </html>
                """
                
                with sync_playwright() as p:
                    browser = p.chromium.launch()
                    page = browser.new_page()
                    page.set_content(html_content)
                    page.wait_for_selector('.mermaid svg', timeout=10000)
                    page.screenshot(path=output_path, full_page=True)
                    browser.close()
                
                return f"Mermaid diagram rendered and saved to: {output_path}"
            except ImportError:
                return f"Error: Could not render mermaid diagram. Playwright not installed. Error details: {str(e)}"
            except Exception as e2:
                return f"Error rendering mermaid diagram: {str(e2)}"
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_initialized = False