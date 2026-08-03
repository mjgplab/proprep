#!/usr/bin/env python3
"""
Transformer Creator - Compatibility Wrapper

DEPRECATED / LEGACY: this is the older "guided v3" transformer creator. It has
been superseded by the interactive table-driven editor
(transformation.table_transformer_creator.TableTransformerCreator), which is the
only transformer-creation path wired into the Redox Site Preparer menu. This
module is retained only for reference and is no longer reachable from the UI; do
not add new callers -- extend TableTransformerCreator instead.

Routes to the v3 guided creator (role labeling, relationship phrases,
demonstrated PDB edits) when a RedoxSite is available, or falls back
to the v2 wizard-based creator otherwise.
"""

import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class TransformerCreator:
    """
    Wrapper class that maintains compatibility with ProPrep's menu system
    and routes to the appropriate creator implementation.

    DEPRECATED / LEGACY: no longer wired into the Redox Site Preparer menu;
    superseded by TableTransformerCreator (see module docstring). Kept for
    reference only.

    When a RedoxSite is provided, uses the v3 guided creator.
    When no RedoxSite is available, falls back to the v2 wizard.
    """

    def __init__(self, output_dir: Optional[Path] = None, processor=None,
                 redox_site=None):
        """
        Initialize transformer creator.

        Args:
            output_dir: Directory to save generated transformers
            processor: ProPrep processor (provides console, session_manager, etc.)
            redox_site: A detected RedoxSite object. If provided, uses v3 creator.
        """
        self.processor = processor
        self.console = processor.console if processor and hasattr(processor, 'console') else None
        self.redox_site = redox_site

        # Determine which creator to use
        if redox_site is not None:
            from .transformer_creator_v3 import TransformerCreatorV3

            self._creator = TransformerCreatorV3(
                redox_site=redox_site,
                output_dir=output_dir,
                console=self.console,
            )
            self._version = "v3"
        else:
            from .transformer_creator_package.creator import TransformerCreator as V2Creator

            # Extract workspace for v2
            workspace = None
            if processor and hasattr(processor, '_get_workspace'):
                try:
                    workspace = processor._get_workspace()
                except Exception as e:
                    logger.warning(f"Could not access workspace: {e}")

            self._creator = V2Creator(output_dir=output_dir, workspace=workspace)
            self._version = "v2"

    def create_transformer(self):
        """
        Run the transformer creation flow.

        Returns True if a transformer was generated, False otherwise.
        """
        try:
            output_file = self._creator.create_transformer()

            if output_file:
                if self.console:
                    self.console.print(f"\n[green]✓ Transformer created successfully![/green]")
                    self.console.print(f"[grey50]File: {output_file}[/grey50]")
                else:
                    print(f"\n✓ Transformer created: {output_file}")
                return True
            else:
                if self.console:
                    self.console.print("\n[yellow]Transformer creation cancelled.[/yellow]")
                else:
                    print("\nTransformer creation cancelled.")
                return False

        except KeyboardInterrupt:
            if self.console:
                self.console.print("\n[yellow]Interrupted by user.[/yellow]")
            else:
                print("\nInterrupted by user.")
            return False

        except Exception as e:
            logger.error(f"Error in transformer creator ({self._version}): {e}")
            if self.console:
                self.console.print(f"\n[red]Error: {e}[/red]")
            else:
                print(f"\nError: {e}")
            import traceback
            traceback.print_exc()
            return False


# For backward compatibility, also allow direct execution
def main():
    """Command-line entry point"""
    from .transformer_creator_package.creator import main as new_main
    new_main()


if __name__ == '__main__':
    main()
