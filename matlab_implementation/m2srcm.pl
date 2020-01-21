#!/usr/bin/perl
if ($#ARGV != 0)
{
  die "Argument must contain filename $#ARGV"
}
else
{
  $fname=$ARGV[0];
}

# If we have a .m file inside a (@)-folder with the same name :
# we will read each file of this folder
if ($fname =~ /^(.*)\@([\d\w-_]*)[\/\\](\2)\.m/)
{
  $name = $2;
  $nameExt = $name.".m";
  $dir = $1."@".$name."/\*.m";
  @fic = glob($dir);
  $i = 0;
  @listeFic[0] = $fname;
  foreach $my_test (@fic)
  {
    if (!($my_test =~ $nameExt))
    {
      $i++;
      @listeFic[$i] = $my_test;
    }
  }
}
# otherwise @-folder, but .m with a different name : ignore it
elsif ($fname =~ /^(.*)\@([\d\w-_]*)[\/\\](.*)\.m/)
{
}
# otherwise
else
{
  @listeFic[0] = $fname;
}
$output = "";
foreach $my_fic (@listeFic)
{

  open(my $in, $my_fic);

  $declTypeDef="";
  $inClass = 0;
  $inAbstractMethodBlock = 0;
  $listeProperties = 0;
  $listeEnumeration = 0;

  $methodAttribute = "";


  while (<$in>)
  {
    if (/(^\s*)(%>)(.*)/)
    {
    }
    elsif(/(^.*)/)
    {
      $output=$output."$1"."\n";
    }
  }
  close $in;
}
print $output;
