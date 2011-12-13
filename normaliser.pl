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
use File::Basename;

#CPAN modules
use Statistics::R;
use Data::Dumper;
use Storable qw(dclone);

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright
printAtStart();

my $options = checkParams();

######################################################################
# CODE HERE
######################################################################
print "Checking if everything makes sense...\n";

# default format is for qiime
my @global_input_format = ("\t",'S','D','','D','M','S','','PN');
if(exists $options->{'if'})
{
    print "Loadind data descriptor from: ". $options->{'if'}. "\n";
    open my $fh, "<", $options->{'if'} or die "**ERROR: Cannot open file: ".$options->{'if'}." $!\n";
    @global_input_format = ();
    while(<$fh>)
    {
        chomp $_;
        if($_ eq ',')
        {
            push @global_input_format, ',';
        }
        else
        {
	        my @fields = split /,/, $_;
	        if($#fields == -1)
	        {
	            push @global_input_format, "";
	        }
	        else
	        {
		        foreach my $char (@fields)
		        {
		            push @global_input_format, $char;
		        }
	        }
        }
    }
    close $fh;
}
my $global_sep = $global_input_format[0];
my $global_type = $global_input_format[-1];
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
$global_data_ref = \@first_data_row;

# we need a place to store the rarefied data
my @global_rarefied_data = ();

# we will need to store the average table too!
my @global_rarefied_average = ();

# an array to hold the distances from the average
my @global_distances = ();

# the index in the rarefied table of the closest to the average
my $global_min_rare_index = 0;

# split the input filename
my($global_tablename, $global_wd, $suffix) = fileparse($options->{'table'});

# get a working directory
if(exists $options->{'working'}) {$global_wd = $options->{'working'}; }
$global_wd .= "/";
`mkdir -p $global_wd`;

# output files
my $global_dist_fh = undef;
if(exists $options->{'dist'}) { open $global_dist_fh, ">", $global_wd.$options->{'dist'} or die "**ERROR: Could not open file: ".$global_wd.$options->{'dist'}." $!\n"; }

my $log_fn = $global_wd.$global_tablename.".log";
if(exists $options->{'log'}) { $log_fn = $global_wd.$options->{'log'}; } 
open my $global_log_fh, ">", $log_fn or die "**ERROR: Could not open file: $log_fn $!\n";

my $global_out_fn_prefix = $global_wd.$global_tablename;
if(exists $options->{'out'}) { $global_out_fn_prefix = $global_wd.$options->{'out'}; } 

# work out the number of reps we need to do
my $global_norm_num_reps = 100;
if(exists $options->{'reps'}) { $global_norm_num_reps = $options->{'reps'}; }

# now, load the DATA table into one huge matrix.
my @global_DATA = (); 
loadDATATable($options->{'table'});

# get the true dimensions of the data
my $global_table_width = $#global_data_matrix;
my $global_table_length = $#{$global_data_matrix[0]};

# make rarefactions
print "Making $global_norm_num_reps multiple rarefactions to $options->{'norm'} individuals...\n";
makeRarefactions();

print "Finding centroid table...\n";
findCentroidTable();

print "Printing output table(s)...\n";
printOutput();

print "Performing Statistical tests...\n";
doStats();

# clean up
print "Cleaning up...\n";
if(exists $options->{'dist'}) {close $global_dist_fh; }
close $global_log_fh;

print "Done!\n----------------------------------------------------------------\n";

######################################################################
# CUSTOM SUBS
######################################################################

sub printOutput
{
    #-----
    # Print the centroid table
    # 
    if(!exists $options->{'of'})
    {
        printRaw();
    }
    else
    {
        my @outfmts = split /,/, $options->{'of'};
        foreach my $fmt (@outfmts)
        {
            if($fmt eq "raw") { printRaw(); }
            elsif($fmt eq "rel") { printRel(); }
            elsif($fmt eq "hel") { printHel(); }
            elsif($fmt eq "bin") { printBin(); }
        }
    }
}

sub printHeader
{
    #-----
    # print the header to an output file
    #
    my ($fh_ref) = @_;
    my $out_fh = ${$fh_ref};
    my $format_index = 1;
    my $sec_index = 0;
    while($global_input_format[$format_index] ne "")
    {
        if($global_input_format[$format_index] eq "D")
        {
            # print the primary descriptor
            print $out_fh join "$global_sep", @primary_column_descriptor;
            print $out_fh "\n";
        }
        elsif($global_input_format[$format_index] eq "S")
        {
            # print a secondary descriptor
            print $out_fh $secondary_column_descriptors[$sec_index]."\n";
            $sec_index++;
        }
        $format_index++;
    }    
    return $format_index;
}

sub printRaw
{
    #-----
    # Print the centroid in raw format (raw integers)
    #
    print "\twriting RAW format\n";
    my $out_fn = $global_out_fn_prefix.".raw.normalised";
    open my $out_fh, ">", $out_fn or die "**ERROR: Could not open file: $out_fn\n";
    
    # get the centroid
    my @centroid = @{$global_rarefied_data[$global_min_rare_index]};

    # first print the headers!
    my $saved_index = printHeader(\$out_fh);
    
    # now print the data
    # we need to keep printing until all the rows are done
    my $row_index = 0;
    $saved_index++;
    foreach my $prd (@primary_row_descriptor)
    {
        my $row_buffer = "";
        my $sec_index = 0;
        foreach my $counter ($saved_index .. $#global_input_format)
        {
            # we need to know the number of info columns
            if($global_input_format[$counter] eq "D")
            {
                # primary descriptor!
                $row_buffer .= $primary_row_descriptor[$row_index]. $global_sep;
            }
            elsif($global_input_format[$counter] eq "S")
            {
                # secondary descriptor!
                $row_buffer .= ${$secondary_row_descriptors[$sec_index]}[$row_index]. $global_sep;
                $sec_index++;
            }
            elsif($global_input_format[$counter] eq "M")
            {
                # data!
                my $cell_index = 0;

                # we either need to print an entire row of the main data
                # matrix or else we need to print one entry from each column
                if($global_type eq 'NP')
                {
                    # print an entire row
                    $row_buffer .= (join $global_sep, @{$centroid[$row_index]}). $global_sep;
                }
                else
                {
                    # $global_type eq 'PN'
                    # print one entry from each column
                    foreach my $i (0 .. $global_table_width)
                    {
                        $row_buffer .=  ${$centroid[$i]}[$row_index].$global_sep;
                    }
                }
            }
            elsif($global_input_format[$counter] eq "")
            {
                last;
            }
        }
        my $row = substr $row_buffer, 0, (length($row_buffer) - 1);
        print $out_fh "$row\n";
        $row_index++;
    }
    close $out_fh;
}

sub printHel
{
    #-----
    # Print the centroid in hellinger format (sqrt of total)
    #
    print "\twriting HEL format\n";
    my $out_fn = $global_out_fn_prefix.".hel.normalised";
    open my $out_fh, ">", $out_fn or die "**ERROR: Could not open file: $out_fn\n";

    # get the centroid
    my @centroid = @{$global_rarefied_data[$global_min_rare_index]};

    # first print the headers!
    my $saved_index = printHeader(\$out_fh);
    
    # now print the data
    # we need to keep printing until all the rows are done
    my $row_index = 0;
    $saved_index++;
    foreach my $prd (@primary_row_descriptor)
    {
        my $row_buffer = "";
        my $sec_index = 0;
        foreach my $counter ($saved_index .. $#global_input_format)
        {
            # we need to know the number of info columns
            if($global_input_format[$counter] eq "D")
            {
                # primary descriptor!
                $row_buffer .= $primary_row_descriptor[$row_index]. $global_sep;
            }
            elsif($global_input_format[$counter] eq "S")
            {
                # secondary descriptor!
                $row_buffer .= ${$secondary_row_descriptors[$sec_index]}[$row_index]. $global_sep;
                $sec_index++;
            }
            elsif($global_input_format[$counter] eq "M")
            {
                # data!
                my $cell_index = 0;

                # we either need to print an entire row of the main data
                # matrix or else we need to print one entry from each column
                if($global_type eq 'NP')
                {
                    # print an entire row
                    foreach my $cell (@{$centroid[$row_index]})
                    {
                        $row_buffer .= sprintf("%.4f",(sqrt $cell)).$global_sep;
                    }
                }
                else
                {
                    # $global_type eq 'PN'
                    # print one entry from each column
                    foreach my $i (0 .. $global_table_width)
                    {
                        $row_buffer .= sprintf("%.4f",sqrt(${$centroid[$i]}[$row_index])).$global_sep;
                    }
                }
            }
            elsif($global_input_format[$counter] eq "")
            {
                last;
            }
        }
        my $row = substr $row_buffer, 0, (length($row_buffer) - 1);
        print $out_fh "$row\n";
        $row_index++;
    }

    close $out_fh;
}

sub printRel
{
    #-----
    # Print the centroid in relative format (percentage of total)
    #
    print "\twriting REL format\n";    
    my $out_fn = $global_out_fn_prefix.".rel.normalised";
    open my $out_fh, ">", $out_fn or die "**ERROR: Could not open file: $out_fn\n";

    # get the centroid
    my @centroid = @{$global_rarefied_data[$global_min_rare_index]};

    # first print the headers!
    my $saved_index = printHeader(\$out_fh);
    
    # now print the data
    # we need to keep printing until all the rows are done
    my $row_index = 0;
    $saved_index++;
    foreach my $prd (@primary_row_descriptor)
    {
        my $row_buffer = "";
        my $sec_index = 0;
        foreach my $counter ($saved_index .. $#global_input_format)
        {
            # we need to know the number of info columns
            if($global_input_format[$counter] eq "D")
            {
                # primary descriptor!
                $row_buffer .= $primary_row_descriptor[$row_index]. $global_sep;
            }
            elsif($global_input_format[$counter] eq "S")
            {
                # secondary descriptor!
                $row_buffer .= ${$secondary_row_descriptors[$sec_index]}[$row_index]. $global_sep;
                $sec_index++;
            }
            elsif($global_input_format[$counter] eq "M")
            {
                # data!
                my $cell_index = 0;

                # we either need to print an entire row of the main data
                # matrix or else we need to print one entry from each column
                if($global_type eq 'NP')
                {
                    # print an entire row
                    foreach my $cell (@{$centroid[$row_index]})
                    {
                        $row_buffer .= sprintf("%.4f",($cell/$options->{'norm'})).$global_sep;
                    }
                }
                else
                {
                    # $global_type eq 'PN'
                    # print one entry from each column
                    foreach my $i (0 .. $global_table_width)
                    {
                        $row_buffer .=  sprintf("%.4f",(${$centroid[$i]}[$row_index]/$options->{'norm'})).$global_sep;
                    }
                }
            }
            elsif($global_input_format[$counter] eq "")
            {
                last;
            }
        }
        my $row = substr $row_buffer, 0, (length($row_buffer) - 1);
        print $out_fh "$row\n";
        $row_index++;
    }
    close $out_fh;
}

sub printBin
{
    #-----
    # Print the centroid in binary format (presence/absense)
    #
    print "\twriting BIN format\n";
    my $out_fn = $global_out_fn_prefix.".bin.normalised";
    open my $out_fh, ">", $out_fn or die "**ERROR: Could not open file: $out_fn\n";

    # get the centroid
    my @centroid = @{$global_rarefied_data[$global_min_rare_index]};

    # first print the headers!
    my $saved_index = printHeader(\$out_fh);
    
    # now print the data
    # we need to keep printing until all the rows are done
    my $row_index = 0;
    $saved_index++;
    foreach my $prd (@primary_row_descriptor)
    {
        my $row_buffer = "";
        my $sec_index = 0;
        foreach my $counter ($saved_index .. $#global_input_format)
        {
            # we need to know the number of info columns
            if($global_input_format[$counter] eq "D")
            {
                # primary descriptor!
                $row_buffer .= $primary_row_descriptor[$row_index]. $global_sep;
            }
            elsif($global_input_format[$counter] eq "S")
            {
                # secondary descriptor!
                $row_buffer .= ${$secondary_row_descriptors[$sec_index]}[$row_index]. $global_sep;
                $sec_index++;
            }
            elsif($global_input_format[$counter] eq "M")
            {
                # data!
                my $cell_index = 0;

                # we either need to print an entire row of the main data
                # matrix or else we need to print one entry from each column
                if($global_type eq 'NP')
                {
                    # print an entire row
                    foreach my $cell (@{$centroid[$row_index]})
                    {
                        if(0 == $cell) { $row_buffer .= "0$global_sep"; }
                        else { $row_buffer .= "1$global_sep"; }
                    }
                }
                else
                {
                    # $global_type eq 'PN'
                    # print one entry from each column
                    foreach my $i (0 .. $global_table_width)
                    {
                        if(0 == ${$centroid[$i]}[$row_index]) { $row_buffer .= "0$global_sep"; }
                        else { $row_buffer .= "1$global_sep"; }
                    }
                }
            }
            elsif($global_input_format[$counter] eq "")
            {
                last;
            }
        }
        my $row = substr $row_buffer, 0, (length($row_buffer) - 1);
        print $out_fh "$row\n";
        $row_index++;
    }
    close $out_fh;
}

sub doStats
{
    #-----
    # Wrapper for doing stats
    #
    # first we print the average and the centroid to two files
    my $ave_fn = $global_wd."AVE.dat";
    my $cent_fn = $global_wd."CENT.dat";
    open my $ave_fh, ">", $ave_fn or die "**ERROR: Cannot open file: $ave_fn for writing $!\n";
    open my $cent_fh, ">", $cent_fn or die "**ERROR: Cannot open file: $cent_fn for writing $!\n";
    
    # get the centroid    
    my @centroid = @{$global_rarefied_data[$global_min_rare_index]};
    
    # and print
    for my $i (0..$global_table_width)
    {
        print $ave_fh ${$global_rarefied_average[$i]}[0];
        print $cent_fh $centroid[$i][0];
        for my $j (1..$global_table_length)
        {
           print $ave_fh ",".${$global_rarefied_average[$i]}[$j];
           print $cent_fh ",".$centroid[$i][$j];
        }
        print $ave_fh "\n";
        print $cent_fh "\n";
    }

    close $ave_fh;
    close $cent_fh;

    if(2 < $global_num_samples)
    {
        #### Create a communication bridge with R and start R
        my $R_instance = Statistics::R->new();
        $R_instance->start();

        # load libraries
        $R_instance->run(qq`library(permute);`);
        $R_instance->run(qq`library(vegan);`);
    
        # read in the tmp files
        my $r_str = "ave<-read.csv(\"$ave_fn\",sep=\",\")";
        $R_instance->run($r_str);  
        print "$r_str\n";  
        $r_str = "cent<-read.csv(\"$cent_fn\",sep=\",\")";
        $R_instance->run($r_str);    
        print "$r_str\n";  
        $R_instance->run(qq`ave<-as.matrix(dist(t(ave), upper=TRUE, diag=TRUE))`);    
        $R_instance->run(qq`cent<-as.matrix(dist(t(cent), upper=TRUE, diag=TRUE))`);
        
        $R_instance->run(qq`mantel.DATA <- mantel(ave,cent);`);
        $R_instance->run(qq`m_stat <- mantel.DATA\$statistic;`);
        $R_instance->run(qq`m_sig <- mantel.DATA\$signif;`);
        print $global_log_fh "Mantel P stat:\t".$R_instance->get('m_sig')."\n";
        print $global_log_fh "Mantel R stat:\t".$R_instance->get('m_stat')."\n";
    }
    else
    {
        print $global_log_fh "Too few samples to perform a statistical tests.\n";
    }
}


sub findCentroidTable
{
    #-----
    # find the rarefied table which is closest to the "average"
    # stored in @global_rarefied_average
    #
    
    # work out the distances from the average
    my %distances = ();
    my $rare_index = 0;
    foreach my $rarefaction (@global_rarefied_data)
    {
        my $dist = findEucDistFromAve($rarefaction);
        $distances{$rare_index} = $dist;
        if(exists $options->{'dist'})
        {
            print $global_dist_fh "$dist\n";
        }
        $rare_index++;
    }

    # find the guy with the least distance, and other stats
    my $max_dist = 0;
    my $min_dist = 1000000000;
    my $mean_dist = 0;
    foreach my $key (keys %distances)
    {
        $mean_dist += $distances{$key};
        if($distances{$key} > $max_dist)
        {
            $max_dist = $distances{$key};
        }

        if($distances{$key} < $min_dist)
        {
            $min_dist = $distances{$key};
            $global_min_rare_index = $key;
        }
    }
    $mean_dist /= $global_norm_num_reps;
    
    print $global_log_fh "---------------------------------------------------\n";
    print $global_log_fh "  Centroid DATA table based normalised statistics\n";
    print $global_log_fh "---------------------------------------------------\n";
    print $global_log_fh "Max dist:\t$max_dist\n";
    print $global_log_fh "Min dist:\t$min_dist\n";
    print $global_log_fh "Range:\t".($max_dist - $min_dist)."\n";
    print $global_log_fh "Mean:\t$mean_dist\n";
}

sub findEucDistFromAve
{
    #-----
    # find the euclidean distance 
    #
    my ($rare_ref) = @_;
    my $dist = 0;
    for my $i (0..$global_table_width)
    {
        for my $j (0..$global_table_length)
        {
            $dist += (${${$rare_ref}[$i]}[$j] - ${$global_rarefied_average[$i]}[$j])**2;
        }
    }
    return sqrt($dist);
}

sub makeRarefactions
{
    #-----
    # make rarefactions of the data table
    # output is an array of rarefied data tables
    # and an "average" table
    #
    # set up our averages martix
    for my $i (0..$global_table_width)
    {
        my @tmp_array = ();
        for my $j (0..$global_table_length)
        {
            push @tmp_array, 0;
        }
        push @global_rarefied_average, \@tmp_array; 
    }
    
    for my $counter (1 .. $global_norm_num_reps)
    {
        my @tmp_data_table = @{ dclone(\@global_data_matrix) };
        rarefyTable(\@tmp_data_table);
        push @global_rarefied_data, \@tmp_data_table;
    }
    
    # now average
    for my $i (0..$global_table_width)
    {
        for my $j (0..$global_table_length)
        {
            ${$global_rarefied_average[$i]}[$j] /= $global_norm_num_reps;
        }
    }
}

sub rarefyTable
{
    #-----
    # Do an individual rarefaction
    # This could be easily parallelised!
    #
    #
    my ($rare_ref) = @_;
    my $row_index = 0;
    foreach my $row (@{$rare_ref})
    {
        # get the total number of individuals
        my $total_individuals = 0;
        
        # use these vars to select which cels to diminish (randomly)
        my %non_zero_indicies = ();
        my $meta_index = 0;
        my $cell_index = 0;
        
        # set up our vars
        foreach my $cell (@{$row})
        {
            if(0 != $cell)
            {
                $total_individuals += $cell;
                $non_zero_indicies{$meta_index} = $cell_index;
                $meta_index++;
            }
            $cell_index++;
        }
        
        # we need to make sure that there are enough guys to make this
        # normalisation feasible
        if($total_individuals < $options->{'norm'})
        {
            if($global_type eq 'PN')
            {
                die "**ERROR: SITE: \"$primary_column_descriptor[$row_index]\" has $total_individuals total entries, and you are trying to normalise to $options->{'norm'}\n";
            }
            else
            {
                die "**ERROR: SITE: \"$primary_row_descriptor[$row_index]\" has $total_individuals total entries, and you are trying to normalise to $options->{'norm'}\n";
            }
        }
        elsif($total_individuals == $options->{'norm'})
        {
            if($global_type eq 'PN')
            {
                print "**WARNING: SITE: \"$primary_column_descriptor[$row_index]\" has $total_individuals total entries, and you are trying to normalise to the same amount\n";
            }
            else
            {
                print "**WARNING: SITE: \"$primary_row_descriptor[$row_index]\" has $total_individuals total entries, and you are trying to normalise to the same amount\n";
            }
        }
        
        # rarefy!
        while($total_individuals > $options->{'norm'})
        {
            my $deduce_index = $non_zero_indicies{int(rand($meta_index))};
            if(0 != ${$row}[$deduce_index])
            {
                ${$row}[$deduce_index]--;
                $total_individuals--;
            }
        }
        
        # add this guy to the averages table
        $cell_index = 0;
        foreach my $cell (@{$row})
        {
            ${$global_rarefied_average[$row_index]}[$cell_index] += $cell;
            $cell_index++;
        }
        $row_index++;
    }
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
        if(1 == $in_header and $global_input_format[$format_index] ne "")
        {
            if($global_input_format[$format_index] eq "D")
            {
                # primary descriptor
                @primary_column_descriptor = split $global_sep, $_;
                $global_total_columns = $#primary_column_descriptor + 1;
                print "Found $global_total_columns columns in total\n";
            }
            elsif($global_input_format[$format_index] eq "S")
            {
                # secondary descriptor 
                push @secondary_column_descriptors, $_;
            }
            else
            {
                die "**ERROR: Unknown character \"".$global_input_format[$format_index]."\" in format file. Exiting.\n";
            }
            $global_num_header_lines++;
        }
        else
        {
            if(1 == $in_header and $global_input_format[$format_index] eq "")
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
                # we are effectively tranposing the array during parsing
                foreach my $field (@data_fields)
                {
                    push @{$global_column_info_indicies[$row_index]}, $field;
                    $row_index++;
                }
            }
        }
    }
    close $fh; 
}

sub createParser
{
    #-----
    # Called on the first data line 
    # sets up the parsing array
    #
    my ($saved_index) = @_;
    my $M_seen = 0;
    foreach my $counter ($saved_index .. $#global_input_format)
    {
        # we need to know the number of info columns
        if($global_input_format[$counter] eq "D" || $global_input_format[$counter] eq "S")
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
        elsif($global_input_format[$counter] eq "M")
        {
            $M_seen = 1;
        }
        elsif($global_input_format[$counter] eq "")
        {
            last;
        }
        else
        {
            die "**ERROR: Unknown character \"".$global_input_format[$counter]."\" in format file. Exiting.\n";
        }
    }
    $global_num_samples = $global_total_columns - $global_num_front_info_columns - $global_num_rear_info_columns - $global_num_ignore_columns;
    print "Data file contains $global_num_front_info_columns front info column(s), $global_num_rear_info_columns rear info column(s) and $global_num_samples data columns\n";
    print "Ignoring $global_num_ignore_columns columns\n----------------------------------------------------------------\n";
    
    # determine the indicies of the data and info columns
    my $row_index = 0;
    foreach my $counter ($saved_index .. $#global_input_format)
    {
        # we need to know the number of info columns
        if($global_input_format[$counter] eq "D")
        {
            # primary descriptor!
            $global_column_info_indicies[$row_index] = \@primary_row_descriptor;
            push @global_is_data_array, 0;
            $row_index++;
        }
        elsif($global_input_format[$counter] eq "S")
        {
            # secondary descriptor!
            my @tmp_desc_array = ();
            $global_column_info_indicies[$row_index] = \@tmp_desc_array;
            push @secondary_row_descriptors, \@tmp_desc_array;
            push @global_is_data_array, 0;
            $row_index++;
        }
        elsif($global_input_format[$counter] eq "M")
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
                $#global_data_matrix = -1; # we need tyo reset this mo-fo
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
        elsif($global_input_format[$counter] eq "")
        {
            last;
        }
        else
        {
            die "**ERROR: Unknown character \"".$global_input_format[$counter]."\" in format file. Exiting.\n";
        }
    }
}

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "of|f:s", "out|o:s", "table|t:s", "norm|n:i", "log|l:s", "dist|d:s", "reps|r:i", "working|w:s", "if|i:s");
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

    # check that the output formats make sense
    if(exists $options{'of'}) 
    {
        my @outfmts = split /,/, $options{'of'};
        foreach my $fmt (@outfmts)
        {
            if(($fmt ne "raw") and ($fmt ne "rel") and ($fmt ne "hel") and ($fmt ne "bin"))
            {
                print "**ERROR: Unknown output format \"$fmt\"\n"; exec("pod2usage $0"); 
            }
        }
    }

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

      -table -t DATA_TABLE                   DATA table to normalise
      -norm -n NORMALISATION_SIZE            Number of sequences to normalise to
      [-reps -r NUM_REPS]                    Number of reps to take (default: 100)
      [-working -w WORKING_DIR]              Place to put any files which will be created (default: location of DATA_TABLE) 
      [-out -o OUT_FILE_PREFIX]              Output normalised file prefix (default: DATA_TABLE)
      [-of -f OUTPUT_TYPE[,OUTPUT_TYPE,...]] Output formats: raw,rel,hel,bin (default: raw only) 
      [-if -i INPUT_FORMAT_FILE]             File to specify the type of the input format (default: QIIME style)                    
      [-log -l LOG_FILE]                     File to store results of mantel tests etc... (default: DATA_TABLE.log)
      [-dist -d DIST_FILE]                   File to store DATA table distances
      [-help -h]                             Displays basic usage information

=cut

