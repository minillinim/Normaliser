#!/usr/bin/perl
###############################################################################
#
#    normaliser.pl
#    
#    Find a centroid DATA table with a given number of individuals
#
#    Copyright (C) 2011 Michael Imelfort and Paul Dennis
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;

#CPAN modules
use Statistics::R;
use File::Basename;
use Text::CSV;
use Data::Dumper;
use Storable qw(dclone);

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get iPNut params and print copyright
printAtStart();

my $options = checkParams();

######################################################################
# CODE HERE
######################################################################
print "Checking if everything makes sense...\n";

# default format is for qiime
my @global_iPNut_format = ('\t','I','D','','D','M','S','','PN');
if(exists $options->{'if'})
{
    print "Loadind data descriptor from: ". $options->{'if'}. "\n";
    open my $fh, "<", $options->{'if'} or die "**ERROR: Cannot open file: ".$options->{'if'}." $!\n";
    @global_iPNut_format = ();
    while(<$fh>)
    {
        chomp $_;
        if($_ eq ',')
        {
            push @global_iPNut_format, ',';
        }
        else
        {
	        my @fields = split /,/, $_;
	        if($#fields == -1)
	        {
	            push @global_iPNut_format, "";
	        }
	        else
	        {
		        foreach my $char (@fields)
		        {
		            push @global_iPNut_format, $char;
		        }
	        }
        }
    }
    close $fh;
}
my $global_sep = $global_iPNut_format[0];
my $global_type = $global_iPNut_format[-1];
print "----------------------------------------------------------------\nData file is of type: $global_type\n";

######################################################################
# csv file particulars
######################################################################
my $global_num_header_lines = 0;
my $global_total_columns = 0;
my $global_num_front_info_columns = 0;
my $global_num_rear_info_columns = 0;
my $global_num_ignore_columns = 0;
my $global_num_samples = 0;
# we need to know which columns in each row are data and which are indicies
my @global_column_info_indicies = ();

# there may be multiple descriptors for each row and column.
# keep track of them using these variables
my @primary_column_descriptor = ();
my @secondary_column_descriptors = ();
my @primary_row_descriptor = ();
my @secondary_row_descriptors = ();

# and finally, somewhere to store the data
my @global_data_matrix = ();
my $global_data_ref;

# for NP type files we need to know whether the column we're looking at contains data
my @global_is_data_array = ();
# grab some memory for the first data row
my @first_data_row = ();
push @global_data_matrix, \@first_data_row;
my $global_data_ref = \@first_data_row;

# we need a place to shove all of the ignores
my @global_no_mans_land = ();

# we need a place to store the rarefied data
my @rarefied_data = ();

# get a working directory
my $global_wd = "working";
if(exists $options->{'working'}) {$global_wd = $options->{'working'}; }
$global_wd .= "/";
`mkdir -p $global_wd`;

# output files
my $global_dist_fh = undef;
if(exists $options->{'dist'}) { open $global_dist_fh, ">", $options->{'dist'} or die "**ERROR: Could not open file: ".$options->{'dist'}." $!\n"; }

my $log_fn = $options->{'table'}.".log";
if(exists $options->{'log'}) { $log_fn = $options->{'log'}; } 
open my $global_log_fh, ">", $log_fn or die "**ERROR: Could not open file: $log_fn $!\n";

my $global_out_fn = $options->{'table'}.".normalised";
if(exists $options->{'out'}) { $global_out_fn = $options->{'out'}; } 

# work out the number of reps we need to do
my $global_norm_num_reps = 1000;
if(exists $options->{'reps'}) { $global_norm_num_reps = $options->{'reps'}; }

# now, load the DATA table into one huge matrix.
my @global_DATA = ();
loadDATATable($options->{'table'});

# parse the table according to the iPNut format

# make rarefactions
print "Making $global_norm_num_reps multiple rarefactions to $options->{'norm'} individuals...\n";
makeRarefactions();

die;
my $cmd_str = "multiple_rarefactions_even_depth.py -i ".$options->{'table'}." -o $global_wd -d ".$options->{'norm'}." -n $global_norm_num_reps --lineages_included --k";
`$cmd_str`;

# find centroid
print "Finding centroid table...\n";
my $centroid_index = find_centroid_table();
my $cp_str = "cp $global_wd"."/rarefaction_".$options->{'norm'}."_$centroid_index".".txt $global_out_fn";
`$cp_str`;

# clean up
print "Cleaning up...\n";
if(exists $options->{'dist'}) {close $global_dist_fh; }
close $global_log_fh;

######################################################################
# CUSTOM SUBS
######################################################################
sub makeRarefactions
{
    #-----
    # make rarefactions of the data table
    #
    my @rarefied_data_table = @{ dclone(\@global_data_matrix) }
    
}
sub loadDATATable
{
    #-----
    # load the DATA table into a large array.
    #
    my ($file_name) = @_;
    my $format_index = 0;
    my $in_header = 1;
    open my $fh, "<", $file_name or die "**ERROR: Cannot open file: $file_name $!\n";
    while(<$fh>)
    {
        chomp $_;
        $format_index++;
        # if we're still in the header, we need to do some special stuff!
        if(1 == $in_header and $global_iPNut_format[$format_index] ne "")
        {
            if($global_iPNut_format[$format_index] eq "I")
            {
                # we ignore this line
            }
            elsif($global_iPNut_format[$format_index] eq "D")
            {
                # primary descriptor
                @primary_column_descriptor = split $global_sep, $_;
                $global_total_columns = $#primary_column_descriptor + 1;
                print "Found $global_total_columns columns in total\n";
            }
            elsif($global_iPNut_format[$format_index] eq "D")
            {
                # secondary descriptor 
                my @tmp_array = split $global_sep, $_;
                push @secondary_column_descriptors, \@tmp_array;
            }
            else
            {
                die "**ERROR: Unknown character \"".$global_iPNut_format[$format_index]."\" in format file. Exiting.\n";
            }
            $global_num_header_lines++;
        }
        else
        {
            if(1 == $in_header and $global_iPNut_format[$format_index] eq "")
            {
                # transition
                $in_header = 0;
                $format_index++;
                # set up the parser
                createParser($format_index);
            }
            elsif($global_type eq 'NP')
            {
                # on the first run, we have a data row in memory
                # on subsequent runs we need to make new ones
                my @tmp_data_row = ();
                push @global_data_matrix, \@tmp_data_row;
                $global_data_ref = \@tmp_data_row;
                # now we definitely have an array in memory
            }
            
            my @data_fields = split $global_sep, $_;
            my $row_index = 0;
            
            # now we're definitely on a data line
            if($global_type eq 'NP')
            {
                # we definitely have an array in memory
                foreach my $field (@data_fields)
                {
                    if(0 == $global_is_data_array[$row_index])
                    {
                        # something!
                        push @{$global_column_info_indicies[$row_index]}, $field;
                    }
                    else
                    {
                        # data!
                        push @{$global_data_ref},$field;
                    }
                    $row_index++;
                }                
            }
            else
            {
                # $global_type eq 'PN'
                # we are effectively transversing the array during parsing
                foreach my $field (@data_fields)
                {
                    push @{$global_column_info_indicies[$row_index]}, $field;
                    $row_index++;
                }
            }
        }
    }
    close $fh; 
    print Dumper @global_data_matrix;
}

sub createParser
{
    #-----
    # Called on the first data line 
    # sets up the parsing array
    #
    my ($saved_index) = @_;
    my $M_seen = 0;
    foreach my $counter ($saved_index .. $#global_iPNut_format)
    {
        # we need to know the number of info columns
        if($global_iPNut_format[$counter] eq "D" || $global_iPNut_format[$counter] eq "S")
        {
            if(0 == $M_seen)
            {
                $global_num_front_info_columns++;
            }
            else
            {
                $global_num_rear_info_columns++;
            }
            
        }
        elsif($global_iPNut_format[$counter] eq "M")
        {
            $M_seen = 1;
        }
        elsif($global_iPNut_format[$counter] eq "I")
        {
            $global_num_ignore_columns++;
        }
        elsif($global_iPNut_format[$counter] eq "")
        {
            last;
        }
        else
        {
            die "**ERROR: Unknown character \"".$global_iPNut_format[$counter]."\" in format file. Exiting.\n";
        }
    }
    $global_num_samples = $global_total_columns - $global_num_front_info_columns - $global_num_rear_info_columns - $global_num_ignore_columns;
    print "Data file contains $global_num_front_info_columns front info column(s), $global_num_rear_info_columns rear info column(s) and $global_num_samples data columns\n";
    print "Ignoring $global_num_ignore_columns columns\n----------------------------------------------------------------\n";
    
    # determine the indicies of the data and info columns
    my $row_index = 0;
    foreach my $counter ($saved_index .. $#global_iPNut_format)
    {
        # we need to know the number of info columns
        if($global_iPNut_format[$counter] eq "D")
        {
            # primary descriptor!
            $global_column_info_indicies[$row_index] = \@primary_row_descriptor;
            push @global_is_data_array, 0;
            $row_index++;
        }
        elsif($global_iPNut_format[$counter] eq "S")
        {
            # secondary descriptor!
            my @tmp_desc_array = ();
            $global_column_info_indicies[$row_index] = \@tmp_desc_array;
            push @secondary_row_descriptors, \@tmp_desc_array;
            push @global_is_data_array, 0;
            $row_index++;
        }
        elsif($global_iPNut_format[$counter] eq "I")
        {
            # Ignore!
            $global_column_info_indicies[$row_index] = \@global_no_mans_land;
            push @global_is_data_array, 0;
            $row_index++;
        }
        elsif($global_iPNut_format[$counter] eq "M")
        {
            # data!
            my $last_row_index = $row_index + $global_num_samples;
            
            # we need to do things a little differently depending on the format of
            # the data. Becuase the dat to normalise may be row-wise or column-wise
            if($global_type eq 'NP')
            {
                # we will normalise row-wise data
	            while($row_index < $last_row_index)
	            {
	                # make all these guys point at a reference to a reference
	                $global_column_info_indicies[$row_index] = \$global_data_ref;
                    push @global_is_data_array, 1;
	                $row_index++;
	            }
            }
            else
            {
                # $global_type eq 'PN'
                while($row_index < $last_row_index)
                {
                    my @tmp_data_column = ();
                    # make all these guys point at a reference to a reference
                    $global_column_info_indicies[$row_index] = \@tmp_data_column;
                    push @global_data_matrix, \@tmp_data_column;
                    $row_index++;
                }
            }
            
        }
        elsif($global_iPNut_format[$counter] eq "")
        {
            last;
        }
        else
        {
            die "**ERROR: Unknown character \"".$global_iPNut_format[$counter]."\" in format file. Exiting.\n";
        }
    }
}

sub find_centroid_table
{
    #-----
    # Find a representative set of $global_norm_sample_size sequences (do this in RRRRRR...)
    #
    #### Create a communication bridge with R and start R
    my $R_instance = Statistics::R->new();
    $R_instance->start();

    my $sampl_p1 = $global_num_samples + 1;
    
    # read in the list of distance matricies
    $R_instance->run(qq`library(foreign);`);
    $R_instance->run(qq`a<-list.files("$global_wd", "*.txt");`);
    
    # work out how many there are and allocate an array
    $R_instance->run(qq`len_a <- length(a);`);
    $R_instance->run(qq`big_frame <- array(0,dim=c($global_num_samples,$global_num_samples,len_a));`);
    
    # load each file individually into a big frame
    my $r_str = "for (i in c(1:len_a)) { j <- i - 1; name <- paste(\"$global_wd\",\"rarefaction_".$options->{'norm'}."\",\"_\",j,\".txt\",sep=\"\"); u<-read.table(name,sep=\"\\t\",row.names=1); u[,$sampl_p1]<-NULL; big_frame[,,i]<-as.matrix(dist(t(u), upper=TRUE, diag=TRUE)); i<-i+1; }";
    $R_instance->run($r_str);    

    # find the average matrix
    $R_instance->run(qq`ave <- big_frame[,,1];`);
    $R_instance->run(qq`for (i in c(2:len_a)) { ave <- ave + big_frame[,,i]; }`);
    $R_instance->run(qq`ave <- ave/len_a;`);
    
    # find the euclidean distance of each matrix from the average
    $R_instance->run(qq`dist<-array(0,dim=c(len_a));`);
    $R_instance->run(qq`for (i in c(1:len_a)) { dist[i] <- sqrt(sum(big_frame[,,i]-ave)^2); }`);
    
    # find the min value
    $R_instance->run(qq`min_index <- which.min(dist);`);
    my $centroid_DATA_index = $R_instance->get('min_index');
    
    # make stats on the distances
    # and log what we did
    $R_instance->run(qq`max_dist <- max(dist);`);
    $R_instance->run(qq`min_dist <- min(dist);`);
    $R_instance->run(qq`range_dist <- max_dist - min_dist;`);
    $R_instance->run(qq`mean_dist <- mean(dist);`);
    $R_instance->run(qq`median_dist <- median(dist);`);
    
    print $global_log_fh "---------------------------------------------------\n";
    print $global_log_fh "  Centroid DATA table based normalised statistics\n";
    print $global_log_fh "---------------------------------------------------\n";
    print $global_log_fh "Max dist:\t".$R_instance->get('max_dist')."\n";
    print $global_log_fh "Min dist:\t".$R_instance->get('min_dist')."\n";
    print $global_log_fh "Range:\t".$R_instance->get('range_dist')."\n";
    print $global_log_fh "Mean:\t".$R_instance->get('mean_dist')."\n";
    print $global_log_fh "Median:\t".$R_instance->get('median_dist')."\n";
  
    if(2 < $global_num_samples)
    {
        $R_instance->run(qq`library(permute);`);
        $R_instance->run(qq`library(vegan);`);
        $R_instance->run(qq`mantel.DATA <- mantel(ave,big_frame[,,min_index]);`);
        $R_instance->run(qq`m_stat <- mantel.DATA\$statistic;`);
        $R_instance->run(qq`m_sig <- mantel.DATA\$signif;`);
        print $global_log_fh "Mantel P stat:\t".$R_instance->get('m_sig')."\n";
        print $global_log_fh "Mantel R stat:\t".$R_instance->get('m_stat')."\n";
    }
    else
    {
        print $global_log_fh "Too few samples to perform a mantel test.\n";
    }
    
    # print all the distances to a file so we can make purdy pictures from them later
    if(exists $options->{'dist'})
    {
        my $num_tables = $R_instance->get('len_a');
        foreach my $counter (1..$num_tables)
        {
            print $global_dist_fh $R_instance->get("dist[$counter]")."\n"         
        }
    }   

    # let the user know the result
    return $centroid_DATA_index - 1;     
}

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "out|o:s", "table|t:s", "norm|n:i", "log|l:s", "dist|d:s", "reps|r:i", "working|w:s", "ot:s", "of:s", "if:s");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsosy items
    #if(!exists $options{''} ) { print "**ERROR: \n"; exec("pod2usage $0"); }
    if(!exists $options{'table'} ) { print "**ERROR: You need to tell me which DATA table to normalise\n"; exec("pod2usage $0"); }
    if(!exists $options{'norm'} ) { print "**ERROR: You need to specify the number of sequences to normalise to\n"; exec("pod2usage $0"); }

    return \%options;
}

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2011 Michael Imelfort and Paul Dennis
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    normaliser.pl
    
=head1 COPYRIGHT

   copyright (C) 2011 Michael Imelfort and Paul Dennis

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

    Find a centroid DATA table with a given number of individuals

=head1 SYNOPSIS

    normaliser.pl -table|t DATA_TABLE -norm|n NORMALISATION_SIZE -log|l LOG_FILE
    
    Normalise a set of DATA tables

      -table -t DATA_TABLE                  DATA table to normalise
      -norm -n NORMALISATION_SIZE           Number of sequences to normalise to
      [-reps -r NUM_REPS]                   Number of reps to take (default: 1000)
      [-working -w WORKING_DIR]             Place to put the multitude of files which will be created (default: working) 
      [-out -o OUT_FILE_PREFIX]             Output normalised file prefix (default: DATA_TABLE.normalised)
      [-ot OUTPUT_TYPE[,OUTPUT_TYPE,...]]   Output types: raw,rel,hel,bin (default: raw only) 
      [-if IPNUT_FORMAT_FILE]               File to specify the type of the iPNut format                    
      [-log -l LOG_FILE]                    File to store results of mantel tests etc... (default: DATA_TABLE.log)
      [-dist -d DIST_FILE]                  File to store DATA table distances
      [-help -h]                            Displays basic usage information

=cut

