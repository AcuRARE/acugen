import click

def_output_dir = '/workspace/acugen/out/'

@click.group()
@click.option('--debug/--no-debug', default=False)
def acugen(debug):
    pass

@acugen.group()

def single():
    pass

@single.command()
@click.option("-i", "--input", "input_file", type=click.Path(file_okay=True, exists=True, resolve_path=True), help="Input file", nargs=1, required=True)
@click.option("-o", "--output", "output_dir", type=click.Path(dir_okay=True, exists=True, resolve_path=True), help="Output directory", nargs=1)
def bam(input_file, output_dir):
    output_dir = def_output_dir if output_dir == None else output_dir
    import acugen_bam
    acugen_bam.main(input_file, output_dir)

@single.command()
@click.option("-i", "--input", "input_file", type=click.Path(file_okay=True, exists=True, resolve_path=True), help="Input files", nargs=2, required=True)
@click.option("-o", "--output", "output_dir", type=click.Path(dir_okay=True, exists=True, resolve_path=True), help="Output directory", nargs=1)
def fastq(input_file, output_dir):
    output_dir = def_output_dir if output_dir == None else output_dir
    import acugen_main
    acugen_main.main(*input_file, output_dir)

@acugen.group()
def batch():
    pass

@batch.command()
@click.option("-i", "--input", "input_dir", type=click.Path(dir_okay=True, exists=True, resolve_path=True), help="Input directory", nargs=1, required=True)
@click.option("-o", "--output", "output_dir", type=click.Path(dir_okay=True, exists=True, resolve_path=True), help="Output directory", nargs=1)
def bam(input_dir, output_dir):
    output_dir = def_output_dir if output_dir == None else output_dir
    import acugen_batch
    acugen_batch.batch_bam(input_dir, output_dir)

@batch.command()
@click.option("-i", "--input", "input_dir", type=click.Path(dir_okay=True, exists=True, resolve_path=True), help="Input directory", nargs=1, required=True)
@click.option("-o", "--output", "output_dir", type=click.Path(dir_okay=True, exists=True, resolve_path=True), help="Output directory", nargs=1)
def fastq(input_dir, output_dir):
    output_dir = def_output_dir if output_dir == None else output_dir
    import acugen_batch
    acugen_batch.batch_fastq(input_dir, output_dir)

if __name__ == '__main__':
    acugen()
