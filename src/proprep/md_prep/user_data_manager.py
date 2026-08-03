"""
User Data Management System

Handles creation, modification, versioning, and metadata tracking for user templates and workflows.
Maintains clear separation between builtin (read-only) and user (writable) content.
"""

import json
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, asdict
import uuid
import hashlib

from rich.console import Console
from rich.table import Table
from rich.panel import Panel

from proprep.utils.paths import get_package_dir
from proprep.utils.prompts import (
    prompt_with_context,
    confirm_with_context,
    int_prompt_with_context,
)


@dataclass
class ContentMetadata:
    """Metadata for user-created or modified content."""
    id: str
    name: str
    description: str
    author: str
    created_date: str
    modified_date: str
    version: str
    content_type: str  # "template" or "workflow"
    source: str  # "custom", "modified_builtin", "builtin"
    parent_id: Optional[str] = None  # For modified content, ID of original
    tags: Optional[List[str]] = None
    
    def __post_init__(self):
        if self.tags is None:
            self.tags = []


@dataclass
class TemplateMetadata(ContentMetadata):
    """Extended metadata for templates."""
    simulation_type: str = ""  # minimization, heating, equilibration, production
    template_path: str = ""
    parameter_count: int = 0
    
    
@dataclass
class WorkflowMetadata(ContentMetadata):
    """Extended metadata for workflows."""
    category: str = ""  # protein, membrane, small_molecule, etc.
    step_count: int = 0
    workflow_path: str = ""
    template_dependencies: Optional[List[str]] = None
    
    def __post_init__(self):
        super().__post_init__()
        if self.template_dependencies is None:
            self.template_dependencies = []


class UserDataManager:
    """Manages user templates, workflows, and metadata with versioning support."""

    def __init__(self, console: Optional[Console] = None):
        self.console = console or Console()

        self.package_dir = get_package_dir()

        # User data always goes to ~/.proprep/ (writable regardless of install method)
        user_data_root = Path.home() / ".proprep"

        # Builtin content is read-only, lives in the package directory
        self.template_base_dir = self.package_dir / "md_templates"
        self.workflow_base_dir = self.package_dir / "md_workflows"

        self.builtin_template_dir = self.template_base_dir / "builtin"
        self.builtin_workflow_dir = self.workflow_base_dir / "builtin"

        # User content is writable, lives in ~/.proprep/
        self.user_template_dir = user_data_root / "md_templates" / "user"
        self.user_workflow_dir = user_data_root / "md_workflows" / "user"

        # Metadata files
        self.template_metadata_file = self.user_template_dir / "metadata.json"
        self.workflow_metadata_file = self.user_workflow_dir / "metadata.json"

        # Ensure directories exist
        self._ensure_directories()

        # Track if we should inform user about data location (for first-time bundle use)
        self._should_show_location_message = False
        if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
            # Only flag to show message if this is a fresh directory (first run)
            if not self.template_metadata_file.exists() and not self.workflow_metadata_file.exists():
                self._should_show_location_message = True

        # Load metadata
        self.template_metadata = self._load_metadata(self.template_metadata_file)
        self.workflow_metadata = self._load_metadata(self.workflow_metadata_file)
        
    def _ensure_directories(self):
        """Ensure all required directories exist."""
        dirs_to_create = [
            self.user_template_dir / "custom",
            self.user_template_dir / "modified",
            self.user_workflow_dir / "custom", 
            self.user_workflow_dir / "modified"
        ]
        
        for directory in dirs_to_create:
            directory.mkdir(parents=True, exist_ok=True)
            
    def _relative_workflow_path(self, workflow_path: Path) -> Path:
        """Compute relative path for a workflow, handling both builtin and user directories."""
        if workflow_path.is_relative_to(self.workflow_base_dir):
            return workflow_path.relative_to(self.workflow_base_dir)
        if workflow_path.is_relative_to(self.user_workflow_dir):
            return Path("user") / workflow_path.relative_to(self.user_workflow_dir)
        return Path(workflow_path.name)

    def _relative_template_path(self, template_path: Path) -> Path:
        """Compute relative path for a template, handling both builtin and user directories."""
        if template_path.is_relative_to(self.template_base_dir):
            return template_path.relative_to(self.template_base_dir)
        if template_path.is_relative_to(self.user_template_dir):
            return Path("user") / template_path.relative_to(self.user_template_dir)
        return Path(template_path.name)

    def _load_metadata(self, metadata_file: Path) -> Dict[str, Dict]:
        """Load metadata from JSON file."""
        if metadata_file.exists():
            try:
                with open(metadata_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                self.console.print(f"[yellow]Warning: Could not load metadata from {metadata_file}: {e}[/yellow]")
                
        return {}
        
    def _save_metadata(self, metadata: Dict[str, Dict], metadata_file: Path):
        """Save metadata to JSON file."""
        try:
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)
        except Exception as e:
            self.console.print(f"[red]Error saving metadata to {metadata_file}: {e}[/red]")
            
    def _generate_content_id(self) -> str:
        """Generate unique content ID."""
        return str(uuid.uuid4())[:8]
        
    def _calculate_content_hash(self, content: str) -> str:
        """Calculate hash of content for change detection."""
        return hashlib.md5(content.encode()).hexdigest()[:8]
        
    def _get_current_timestamp(self) -> str:
        """Get current timestamp string."""
        return datetime.now().isoformat()

    def _show_location_message_if_needed(self):
        """Show user data location message on first use (for bundled executables)."""
        if self._should_show_location_message:
            user_data_root = Path.home() / ".proprep"
            self.console.print(
                f"\n[cyan]User data will be stored in: {user_data_root}[/cyan]",
                style="grey50"
            )
            self._should_show_location_message = False  # Only show once

    def create_template(self, name: str, description: str, simulation_type: str,
                       content: str, author: str = "User", tags: Optional[List[str]] = None) -> str:
        """Create a new user template."""

        # Show location message on first use
        self._show_location_message_if_needed()
        
        # Validate simulation type
        valid_types = ["minimization", "heating", "equilibration", "production"]
        if simulation_type not in valid_types:
            raise ValueError(f"Invalid simulation type. Must be one of: {valid_types}")
            
        # Generate ID and paths
        template_id = self._generate_content_id()
        # Truncate name portion to avoid filesystem limits
        name_for_filename = name.lower().replace(' ', '_')
        if len(name_for_filename) > 60:
            name_for_filename = name_for_filename[:60]
        template_filename = f"{template_id}_{name_for_filename}.mdin"
        template_path = self.user_template_dir / "custom" / template_filename
        
        # Save template content
        with open(template_path, 'w') as f:
            f.write(content)
            
        # Create metadata
        metadata = TemplateMetadata(
            id=template_id,
            name=name,
            description=description,
            author=author,
            created_date=self._get_current_timestamp(),
            modified_date=self._get_current_timestamp(),
            version="1.0",
            content_type="template",
            source="custom",
            tags=tags or [],
            simulation_type=simulation_type,
            template_path=str(self._relative_template_path(template_path)),
            parameter_count=self._count_parameters(content)
        )
        
        # Save metadata
        self.template_metadata[template_id] = asdict(metadata)
        self._save_metadata(self.template_metadata, self.template_metadata_file)
        
        self.console.print(f"[green]✓ Template '{name}' created with ID: {template_id}[/green]")
        return template_id
        
    def modify_template(self, source_path: str, name: str, description: str, 
                       content: str, author: str = "User") -> str:
        """Create a modified version of an existing template."""
        
        # Determine if source is builtin or user template
        source_id = None
        is_builtin = False
        
        # Check if it's a builtin template (by path)
        if source_path.startswith("builtin/"):
            is_builtin = True
            parent_id = f"builtin_{source_path.replace('/', '_').replace('.mdin', '')}"
        else:
            # Find source template ID
            for tid, meta in self.template_metadata.items():
                if meta.get("template_path") == source_path:
                    source_id = tid
                    parent_id = tid
                    break
            else:
                raise ValueError(f"Source template not found: {source_path}")
                
        # Generate new ID and path
        template_id = self._generate_content_id()
        # Truncate name portion to avoid filesystem limits
        name_for_filename = name.lower().replace(' ', '_')
        if len(name_for_filename) > 50:  # Slightly shorter for "modified_" prefix
            name_for_filename = name_for_filename[:50]
        template_filename = f"{template_id}_modified_{name_for_filename}.mdin"
        template_path = self.user_template_dir / "modified" / template_filename
        
        # Save modified content
        with open(template_path, 'w') as f:
            f.write(content)
            
        # Create metadata
        metadata = TemplateMetadata(
            id=template_id,
            name=f"Modified: {name}",
            description=f"Modified version of '{name}': {description}",
            author=author,
            created_date=self._get_current_timestamp(),
            modified_date=self._get_current_timestamp(),
            version="1.0",
            content_type="template",
            source="modified_builtin" if is_builtin else "modified_user",
            parent_id=parent_id,
            simulation_type=self._determine_simulation_type(source_path),
            template_path=str(self._relative_template_path(template_path)),
            parameter_count=self._count_parameters(content)
        )
        
        # Save metadata
        self.template_metadata[template_id] = asdict(metadata)
        self._save_metadata(self.template_metadata, self.template_metadata_file)
        
        self.console.print(f"[green]✓ Modified template saved with ID: {template_id}[/green]")
        return template_id
        
    def create_workflow(self, name: str, description: str, category: str,
                       steps: List[Dict], author: str = "User", tags: List[str] = None) -> str:
        """Create a new user workflow."""

        # Show location message on first use
        self._show_location_message_if_needed()

        # Generate ID and paths
        workflow_id = self._generate_content_id()
        # Truncate name portion to avoid filesystem limits
        name_for_filename = name.lower().replace(' ', '_')
        if len(name_for_filename) > 60:
            name_for_filename = name_for_filename[:60]
        workflow_filename = f"{workflow_id}_{name_for_filename}.json"
        workflow_path = self.user_workflow_dir / "custom" / workflow_filename
        
        # Create workflow data
        workflow_data = {
            "name": name,
            "description": description,
            "version": "1.0",
            "author": author,
            "category": category,
            "steps": steps
        }
        
        # Save workflow
        with open(workflow_path, 'w') as f:
            json.dump(workflow_data, f, indent=2)
            
        # Extract template dependencies
        template_deps = []
        for step in steps:
            if "template" in step:
                template_deps.append(step["template"])
                
        # Create metadata
        metadata = WorkflowMetadata(
            id=workflow_id,
            name=name,
            description=description,
            author=author,
            created_date=self._get_current_timestamp(),
            modified_date=self._get_current_timestamp(),
            version="1.0",
            content_type="workflow",
            source="custom",
            tags=tags or [],
            category=category,
            step_count=len(steps),
            workflow_path=str(self._relative_workflow_path(workflow_path)),
            template_dependencies=template_deps
        )

        # Save metadata
        self.workflow_metadata[workflow_id] = asdict(metadata)
        self._save_metadata(self.workflow_metadata, self.workflow_metadata_file)

        self.console.print(f"[green]✓ Workflow '{name}' created with ID: {workflow_id}[/green]")
        return workflow_id
        
    def modify_workflow(self, source_path: str, name: str, description: str, 
                       steps: List[Dict], author: str = "User") -> str:
        """Create a modified version of an existing workflow."""
        
        # Determine source
        is_builtin = source_path.startswith("builtin/")
        if is_builtin:
            parent_id = f"builtin_{Path(source_path).stem}"
        else:
            # Find source workflow ID
            parent_id = None
            for wid, meta in self.workflow_metadata.items():
                if meta.get("workflow_path") == source_path:
                    parent_id = wid
                    break
            if not parent_id:
                raise ValueError(f"Source workflow not found: {source_path}")
                
        # Generate new ID and path
        workflow_id = self._generate_content_id()
        # Truncate name portion to avoid filesystem limits
        name_for_filename = name.lower().replace(' ', '_')
        if len(name_for_filename) > 50:  # Slightly shorter for "modified_" prefix
            name_for_filename = name_for_filename[:50]
        workflow_filename = f"{workflow_id}_modified_{name_for_filename}.json"
        workflow_path = self.user_workflow_dir / "modified" / workflow_filename
        
        # Create workflow data
        workflow_data = {
            "name": f"Modified: {name}",
            "description": f"Modified version of '{name}': {description}",
            "version": "1.0",
            "author": author,
            "category": "custom",
            "steps": steps
        }
        
        # Save workflow
        with open(workflow_path, 'w') as f:
            json.dump(workflow_data, f, indent=2)
            
        # Extract template dependencies
        template_deps = []
        for step in steps:
            if "template" in step:
                template_deps.append(step["template"])
                
        # Create metadata
        metadata = WorkflowMetadata(
            id=workflow_id,
            name=f"Modified: {name}",
            description=f"Modified version of '{name}': {description}",
            author=author,
            created_date=self._get_current_timestamp(),
            modified_date=self._get_current_timestamp(),
            version="1.0",
            content_type="workflow",
            source="modified_builtin" if is_builtin else "modified_user",
            parent_id=parent_id,
            category="custom",
            step_count=len(steps),
            workflow_path=str(self._relative_workflow_path(workflow_path)),
            template_dependencies=template_deps
        )
        
        # Save metadata
        self.workflow_metadata[workflow_id] = asdict(metadata)
        self._save_metadata(self.workflow_metadata, self.workflow_metadata_file)
        
        self.console.print(f"[green]✓ Modified workflow saved with ID: {workflow_id}[/green]")
        return workflow_id
        
    def list_templates(self, show_builtin: bool = True, show_user: bool = True) -> Dict[str, Dict]:
        """List available templates with metadata."""
        templates = {}
        
        if show_builtin:
            # Add builtin templates (scan protocol directories dynamically)
            for protocol_dir in sorted(self.builtin_template_dir.iterdir()):
                if not protocol_dir.is_dir():
                    continue
                for template_file in sorted(protocol_dir.glob("*.mdin")):
                    sim_type = self._determine_simulation_type(template_file.name)
                    template_id = f"builtin_{protocol_dir.name}_{template_file.stem}"

                    # Try to parse template header for metadata
                    try:
                        with open(template_file, 'r') as f:
                            content = f.read()
                        header_metadata = self._parse_template_header(content)
                    except Exception:
                        header_metadata = {}

                    # Use header metadata if available, otherwise fall back to defaults
                    name = header_metadata.get('name', template_file.stem.replace('_', ' ').title())
                    description = header_metadata.get('description', f"Built-in {sim_type} template")
                    author = header_metadata.get('author', "ProPrep")
                    priority = header_metadata.get('priority', 999)  # Default low priority
                    source_url = header_metadata.get('source_url', None)

                    templates[template_id] = {
                        "id": template_id,
                        "name": name,
                        "description": description,
                        "author": author,
                        "priority": priority,
                        "source": "builtin",
                        "source_url": source_url,
                        "simulation_type": sim_type,
                        "template_path": f"builtin/{protocol_dir.name}/{template_file.name}",
                        "content_type": "template"
                    }
                        
        if show_user:
            # Add user templates from metadata
            for template_id, metadata in self.template_metadata.items():
                templates[template_id] = metadata

            # Also scan custom directory for JSON templates not in metadata
            custom_dir = self.user_template_dir / "custom"
            if custom_dir.exists():
                for template_file in custom_dir.glob("*.json"):
                    try:
                        with open(template_file, 'r') as f:
                            template_data = json.load(f)
                            template_id = template_data.get('id', template_file.stem)
                            # Only add if not already in templates (metadata takes precedence)
                            if template_id not in templates:
                                templates[template_id] = template_data
                    except (json.JSONDecodeError, IOError):
                        continue

        return templates
        
    def list_workflows(self, show_builtin: bool = True, show_user: bool = True) -> Dict[str, Dict]:
        """List available workflows with metadata."""
        workflows = {}
        
        if show_builtin:
            # Add builtin workflows (synthesized metadata)
            for workflow_file in sorted(self.builtin_workflow_dir.glob("*.json")):
                try:
                    with open(workflow_file, 'r') as f:
                        workflow_data = json.load(f)
                        
                    workflow_id = f"builtin_{workflow_file.stem}"
                    workflows[workflow_id] = {
                        "id": workflow_id,
                        "name": workflow_data.get("name", workflow_file.stem),
                        "description": workflow_data.get("description", "Built-in workflow"),
                        "author": workflow_data.get("author", "ProPrep"),
                        "source": "builtin",
                        "category": workflow_data.get("category", "protein"),
                        "step_count": len(workflow_data.get("steps", [])),
                        "workflow_path": f"builtin/{workflow_file.name}",
                        "content_type": "workflow"
                    }
                except Exception as e:
                    self.console.print(f"[yellow]Warning: Could not load builtin workflow {workflow_file}: {e}[/yellow]")
                    
        if show_user:
            # Add user workflows from metadata
            for workflow_id, metadata in self.workflow_metadata.items():
                workflows[workflow_id] = metadata
                
        return workflows
        
    def get_template_content(self, template_id_or_path: str) -> Tuple[str, Dict]:
        """Get template content and metadata."""
        if template_id_or_path.startswith("builtin/"):
            # Handle builtin template by path
            template_path = self.template_base_dir / template_id_or_path
            if not template_path.exists():
                raise FileNotFoundError(f"Template not found: {template_path}")
                
            with open(template_path, 'r') as f:
                content = f.read()
                
            # Synthesize metadata for builtin templates
            sim_type = self._determine_simulation_type(template_id_or_path)

            metadata = {
                "name": template_path.stem.replace('_', ' ').title(),
                "description": f"Built-in {sim_type} template",
                "source": "builtin",
                "simulation_type": sim_type
            }
            
            return content, metadata
            
        elif template_id_or_path in self.template_metadata:
            # Handle user template by ID
            metadata = self.template_metadata[template_id_or_path]
            template_path = self.template_base_dir / metadata["template_path"]
            
            if not template_path.exists():
                raise FileNotFoundError(f"Template file not found: {template_path}")
                
            with open(template_path, 'r') as f:
                content = f.read()
                
            return content, metadata
            
        else:
            raise ValueError(f"Template not found: {template_id_or_path}")
            
    def get_workflow_content(self, workflow_id_or_path: str) -> Tuple[Dict, Dict]:
        """Get workflow content and metadata."""
        if workflow_id_or_path.startswith("builtin/"):
            # Handle builtin workflow by path
            workflow_path = self.workflow_base_dir / workflow_id_or_path
            if not workflow_path.exists():
                raise FileNotFoundError(f"Workflow not found: {workflow_path}")
                
            with open(workflow_path, 'r') as f:
                content = json.load(f)
                
            # Synthesize metadata for builtin workflows
            metadata = {
                "name": content.get("name", workflow_path.stem),
                "description": content.get("description", "Built-in workflow"),
                "source": "builtin",
                "category": content.get("category", "protein")
            }
            
            return content, metadata
            
        elif workflow_id_or_path in self.workflow_metadata:
            # Handle user workflow by ID
            metadata = self.workflow_metadata[workflow_id_or_path]
            workflow_path = self.workflow_base_dir / metadata["workflow_path"]
            
            if not workflow_path.exists():
                raise FileNotFoundError(f"Workflow file not found: {workflow_path}")
                
            with open(workflow_path, 'r') as f:
                content = json.load(f)
                
            return content, metadata
            
        else:
            raise ValueError(f"Workflow not found: {workflow_id_or_path}")
            
    def delete_user_content(self, content_id: str) -> bool:
        """Delete user template or workflow."""
        # Check templates
        if content_id in self.template_metadata:
            metadata = self.template_metadata[content_id]
            template_path = self.template_base_dir / metadata["template_path"]
            
            if template_path.exists():
                template_path.unlink()
                
            del self.template_metadata[content_id]
            self._save_metadata(self.template_metadata, self.template_metadata_file)
            
            self.console.print(f"[green]✓ Template '{metadata['name']}' deleted[/green]")
            return True
            
        # Check workflows
        if content_id in self.workflow_metadata:
            metadata = self.workflow_metadata[content_id]
            workflow_path = self.workflow_base_dir / metadata["workflow_path"]
            
            if workflow_path.exists():
                workflow_path.unlink()
                
            del self.workflow_metadata[content_id]
            self._save_metadata(self.workflow_metadata, self.workflow_metadata_file)
            
            self.console.print(f"[green]✓ Workflow '{metadata['name']}' deleted[/green]")
            return True
            
        self.console.print(f"[yellow]Content not found: {content_id}[/yellow]")
        return False
        
    def export_content(self, content_id: str, export_path: Path) -> bool:
        """Export template or workflow for sharing."""
        try:
            # Check if it's a template
            if content_id in self.template_metadata:
                metadata = self.template_metadata[content_id]
                template_path = self.template_base_dir / metadata["template_path"]
                
                export_data = {
                    "type": "template",
                    "metadata": metadata,
                    "content": template_path.read_text()
                }
                
            # Check if it's a workflow
            elif content_id in self.workflow_metadata:
                metadata = self.workflow_metadata[content_id]
                workflow_path = self.workflow_base_dir / metadata["workflow_path"]
                
                with open(workflow_path, 'r') as f:
                    workflow_content = json.load(f)
                    
                export_data = {
                    "type": "workflow",
                    "metadata": metadata,
                    "content": workflow_content
                }
                
            else:
                self.console.print(f"[red]Content not found: {content_id}[/red]")
                return False
                
            # Save export file
            with open(export_path, 'w') as f:
                json.dump(export_data, f, indent=2)
                
            self.console.print(f"[green]✓ Content exported to {export_path}[/green]")
            return True
            
        except Exception as e:
            self.console.print(f"[red]Export failed: {e}[/red]")
            return False
            
    def import_content(self, import_path: Path, author: str = "User") -> Optional[str]:
        """Import template or workflow from export file."""
        try:
            with open(import_path, 'r') as f:
                import_data = json.load(f)
                
            content_type = import_data.get("type")
            metadata = import_data.get("metadata", {})
            content = import_data.get("content")
            
            if content_type == "template":
                return self.create_template(
                    name=metadata.get("name", "Imported Template"),
                    description=f"Imported: {metadata.get('description', '')}",
                    simulation_type=metadata.get("simulation_type", "minimization"),
                    content=content,
                    author=author,
                    tags=metadata.get("tags", [])
                )
                
            elif content_type == "workflow":
                return self.create_workflow(
                    name=metadata.get("name", "Imported Workflow"),
                    description=f"Imported: {metadata.get('description', '')}",
                    category=metadata.get("category", "custom"),
                    steps=content.get("steps", []),
                    author=author,
                    tags=metadata.get("tags", [])
                )
                
            else:
                self.console.print(f"[red]Unknown content type: {content_type}[/red]")
                return None
                
        except Exception as e:
            self.console.print(f"[red]Import failed: {e}[/red]")
            return None
            
    def _count_parameters(self, content: str) -> int:
        """Count AMBER parameters in content."""
        import re
        # Simple parameter counting - count lines with "param = value,"
        return len(re.findall(r'^\s*\w+\s*=\s*[^,\n]+,', content, re.MULTILINE))
        
    def _parse_template_header(self, content: str) -> Dict[str, str]:
        """Extract metadata from template header comments."""
        metadata = {}
        lines = content.split('\n')
        
        for line in lines[:10]:  # Only check first 10 lines
            line = line.strip()
            if line.startswith('! TEMPLATE:'):
                metadata['name'] = line.replace('! TEMPLATE:', '').strip()
            elif line.startswith('! DESCRIPTION:'):
                metadata['description'] = line.replace('! DESCRIPTION:', '').strip()
            elif line.startswith('! TYPE:'):
                metadata['type'] = line.replace('! TYPE:', '').strip()
            elif line.startswith('! AUTHOR:'):
                metadata['author'] = line.replace('! AUTHOR:', '').strip()
            elif line.startswith('! VERSION:'):
                metadata['version'] = line.replace('! VERSION:', '').strip()
            elif line.startswith('! PRIORITY:'):
                try:
                    metadata['priority'] = int(line.replace('! PRIORITY:', '').strip())
                except ValueError:
                    metadata['priority'] = 999  # Default low priority for invalid values
            elif line.startswith('! SOURCE:'):
                metadata['source_url'] = line.replace('! SOURCE:', '').strip()
            elif line.startswith('&cntrl') or line.startswith('&'):
                break  # Stop at first control section
                
        return metadata
        
    def _determine_simulation_type(self, template_path: str) -> str:
        """Determine simulation type from template path."""
        if "minimization" in template_path:
            return "minimization"
        elif "heating" in template_path:
            return "heating"  
        elif "equilibration" in template_path:
            return "equilibration"
        elif "production" in template_path:
            return "production"
        else:
            return "unknown"
            
    def save_custom_workflow(self, workflow_data: Dict) -> str:
        """Save custom workflow and return its ID."""
        workflow_id = workflow_data.get('id', str(uuid.uuid4()))
        workflow_name = workflow_data.get('name', f'workflow_{workflow_id}')

        # Create filename with truncated name to avoid filesystem limits
        # UUID is 36 chars + underscore + .json (5 chars) = 42 chars overhead
        # Leave room for ~60 char name portion to stay well under 255 char limit
        name_for_filename = workflow_name.lower().replace(' ', '_')
        if len(name_for_filename) > 60:
            name_for_filename = name_for_filename[:60]
        workflow_filename = f"{workflow_id}_{name_for_filename}.json"
        workflow_path = self.user_workflow_dir / "custom" / workflow_filename
        
        # Ensure directory exists
        workflow_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save workflow data
        with open(workflow_path, 'w') as f:
            json.dump(workflow_data, f, indent=2)
        
        # Update metadata dictionary so list_workflows() can find it
        self.workflow_metadata[workflow_id] = {
            "name": workflow_name,
            "description": workflow_data.get('description', ''),
            "source": "custom",
            "category": workflow_data.get('category', 'custom'),
            "created_date": workflow_data.get('created_date', datetime.now().isoformat()),
            "modified_date": workflow_data.get('modified_date', datetime.now().isoformat()),
            "workflow_path": str(self._relative_workflow_path(workflow_path)),
            "content_type": "workflow"
        }
        
        # Save metadata file
        self._save_metadata(self.workflow_metadata, self.workflow_metadata_file)
            
        return workflow_id
        
    def load_custom_workflow(self, workflow_id: str) -> Optional[Dict]:
        """Load custom workflow by ID."""
        # Search in custom workflows
        custom_dir = self.user_workflow_dir / "custom"
        if custom_dir.exists():
            for workflow_file in custom_dir.glob(f"{workflow_id}_*.json"):
                try:
                    with open(workflow_file, 'r') as f:
                        return json.load(f)
                except (json.JSONDecodeError, IOError):
                    continue
        return None
        
    def list_saved_protocols(self) -> Dict[str, Dict]:
        """List user-saved protocols (convenience wrapper over list_workflows)."""
        return self.list_workflows(show_builtin=False, show_user=True)

    def save_custom_template(self, template_data: Dict, template_id: str = None) -> str:
        """Save custom template and return its ID."""
        if template_id is None:
            template_id = str(uuid.uuid4())

        template_name = template_data.get('name', f'template_{template_id}')

        # Create filename with truncated name to avoid filesystem limits
        # UUID is 36 chars + underscore + .json (5 chars) = 42 chars overhead
        # Leave room for ~60 char name portion to stay well under 255 char limit
        name_for_filename = template_name.lower().replace(' ', '_')
        if len(name_for_filename) > 60:
            name_for_filename = name_for_filename[:60]
        template_filename = f"{template_id}_{name_for_filename}.json"
        template_path = self.user_template_dir / "custom" / template_filename
        
        # Ensure directory exists
        template_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Add ID to template data if not present
        template_data['id'] = template_id
        
        # Save template data
        with open(template_path, 'w') as f:
            json.dump(template_data, f, indent=2)
            
        return template_id
        
    def load_custom_template(self, template_id: str) -> Optional[Dict]:
        """Load custom template by ID."""
        # Search in custom templates
        custom_dir = self.user_template_dir / "custom"
        if custom_dir.exists():
            for template_file in custom_dir.glob(f"{template_id}_*.json"):
                try:
                    with open(template_file, 'r') as f:
                        return json.load(f)
                except (json.JSONDecodeError, IOError):
                    continue
        return None

    def get_workflow_template_ids(self, workflow_id: str) -> set:
        """
        Extract all template references from a workflow.

        Checks both new-style template_ref (builtin paths) and legacy
        custom_template_id (UUIDs) fields.

        Args:
            workflow_id: ID of the workflow to analyze

        Returns:
            Set of template IDs/refs used by the workflow
        """
        template_ids = set()

        try:
            workflow_data = self.load_custom_workflow(workflow_id)
            if workflow_data:
                for step in workflow_data.get('steps', []):
                    # New-style: template_ref (builtin path)
                    template_ref = step.get('template_ref')
                    if template_ref:
                        template_ids.add(template_ref)
                    # Legacy: custom_template_id (UUID)
                    custom_template_id = step.get('custom_template_id')
                    if custom_template_id:
                        template_ids.add(custom_template_id)
        except Exception:
            # If workflow can't be loaded, return empty set
            pass

        return template_ids

    def find_template_usage(self, template_id: str) -> list:
        """
        Find all workflows that reference a specific template.

        Args:
            template_id: ID of the template to search for

        Returns:
            List of workflow IDs that use this template
        """
        used_by_workflows = []

        # Get all custom workflows
        workflows = self.list_workflows(show_builtin=False, show_user=True)

        for workflow_id in workflows.keys():
            try:
                workflow_data = self.load_custom_workflow(workflow_id)
                if workflow_data:
                    # Check if any step uses this template (new or legacy field)
                    for step in workflow_data.get('steps', []):
                        if step.get('template_ref') == template_id or step.get('custom_template_id') == template_id:
                            used_by_workflows.append(workflow_id)
                            break  # No need to check other steps
            except Exception:
                # Skip workflows that can't be loaded
                continue

        return used_by_workflows