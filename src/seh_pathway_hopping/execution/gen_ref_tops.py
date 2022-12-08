import click

@click.command()
@click.argument('gexp')
def gen_ref_tops(gexp):


    from seh_pathway_hopping._tasks import (
      write_lig_ref_selection_pdbs,
      write_lig_centered_ref_selection_pdbs,
      GEXP_LIG_IDS
    )

    lig_id = dict(GEXP_LIG_IDS)[gexp]

    write_lig_ref_selection_pdbs(lig_id)
    write_lig_centered_ref_selection_pdbs(lig_id)

if __name__ == "__main__":

    gen_ref_tops()
