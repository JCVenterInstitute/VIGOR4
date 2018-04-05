package org.jcvi.vigor.service;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.vigor.AppConfig;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.VigorUtils;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.CoreMatchers.hasItem;
import static org.junit.Assert.assertThat;

@RunWith(SpringRunner.class)
@ContextConfiguration(classes = AppConfig.class)


public class PeptideServiceTest {

    private final static Logger LOGGER = LogManager.getLogger(PeptideServiceTest.class);

    @Autowired
    private PeptideService peptideService;

    @Test
    public void testPeptides() throws ServiceException {

        String[][] testingData = new String[][]{
                new String[]{
                        "HA_1_615",
                        String.join("",
                                "MKAILVVLLYTFATANADTLCIGYHANNST",
                                "DTVDTVLEKNVTVTHSVNLLEDKHNGKLCK",
                                "LRGVAPLHLGKCNIAGWILGNPECESLSTA",
                                "SSWSYIVETSSSDNGTCYPGDFINYEELRE",
                                "QLSSVSSFERFEIFPKTSSWPNHDSNKGVT",
                                "AACPQAGAKSFYKNLIWLVKKGNSYPKLSK",
                                "SYINDKGKEVLVLWGIHH"),
                        "flua_ha_mp",
                        "MKAILVVLLYTFATANA",
                        String.join("",
                                "DTLCIGYHANNSTDTVDTVLEKNVTVTHSV",
                                "NLLEDKHNGKLCKLRGVAPLHLGKCNIAGW",
                                "ILGNPECESLSTASSWSYIVETSSSDNGTC",
                                "YPGDFINYEELREQLSSVSSFERFEIFPKT",
                                "SSWPNHDSNKGVTAACPQAGAKSFYKNLIW",
                                "LVKKGNSYPKLSKSYINDKGKEVLVLWGIHH")
                },
                new String[] {
                        "HA_1129_1753.1",
                        String.join("",
                                "QNAQGSGYAADLKSTQNAIDKITNKVNSVI",
                                "EKMNTQFTAVGKEFNHLEKRIENLNKKVDD",
                                "GFLDIWTYNAELLVLLENERTLDYHDSNVK",
                                "NLYEKVRSQLKNNAKEIGNGCFEFYHKCDN",
                                "TCMESVKNGTYDYPRYSEEAKLNREEIDGV",
                                "KLESTRIYQILAIYSTVASSLVLIVSLGAI",
                                "SFWMCSNGSLQCRICI*"),
                        "flua_ha_mp",
                        String.join("",
                                "QNAQGSGYAADLKSTQNAIDKITNKVNSVI",
                                "EKMNTQFTAVGKEFNHLEKRIENLNKKVDD",
                                "GFLDIWTYNAELLVLLENERTLDYHDSNVK",
                                "NLYEKVRSQLKNNAKEIGNGCFEFYHKCDN",
                                "TCMESVKNGTYDYPRYSEEAKLNREEIDGV",
                                "KLESTRIYQILAIYSTVASSLVLIVSLGAI",
                                "SFWMCSNGSLQCRICI")
                },

                new String[] {
                        "gi_260907762.1",
                        String.join("",
                                "MNIQILVFALVAIIPTNADKICLGHHAVSN",
                                "GTKVNTLTERGVEVVNATETVERTNVPRIC",
                                "SKGKRTVDLGQCGLLGTITGPPQCDQFLEF",
                                "SADLIIERRGGSDVCYPGKFVNEEALRQIL",
                                "RESGGIDKETMGFTYSGIRTNGATSACRRS",
                                "GSSFYAEMKWLLSNTDNAAFPQMTKSYKNT",
                                "RKDPALIIWGIHHSGSTTEQTKLYGSGSKL",
                                "ITVGSSNYQQSFVPSPGARPQVNGQSGRID",
                                "FHWLMLNPNDTVTFSFNGAFIAPDRASFLK",
                                "GKSMGIQSGVQVDANCEGDCYHSGGTIISN",
                                "LPFQNINSRAVGKCPRYVKQESLMLATGMK",
                                "NVPELPKGRGLFGAIAGFIENGWEGLIDGW",
                                "YGFRHQNAQGEGTAADYKSTQSAIDQITGK",
                                "LNRLIEKTNQQFELIDNEFTEVEKQIGNVI",
                                "NWTRDSMTEVWSYNAELLVAMENQHTIDLA",
                                "DSEMNKLYERVKRQLRENAEEDGTGCFEIF",
                                "HKCDDDCMASIRNNTYDHSKYREEAMQNRI",
                                "QINPVKLSSGYKDVILWFSFGASCFILLAI",
                                "AMGLVFICVKNGNMRCTICI*"),
                        "flua_ha_mp",
                        "MNIQILVFALVAIIPTNA",
                        String.join("",
                                "DKICLGHHAVSNGTKVNTLTERGVEVVNAT",
                                "ETVERTNVPRICSKGKRTVDLGQCGLLGTI",
                                "TGPPQCDQFLEFSADLIIERRGGSDVCYPG",
                                "KFVNEEALRQILRESGGIDKETMGFTYSGI",
                                "RTNGATSACRRSGSSFYAEMKWLLSNTDNA",
                                "AFPQMTKSYKNTRKDPALIIWGIHHSGSTT",
                                "EQTKLYGSGSKLITVGSSNYQQSFVPSPGA",
                                "RPQVNGQSGRIDFHWLMLNPNDTVTFSFNG",
                                "AFIAPDRASFLKGKSMGIQSGVQVDANCEG",
                                "DCYHSGGTIISNLPFQNINSRAVGKCPRYV",
                                "KQESLMLATGMKNVPELPKGR"),
                        String.join("",
                                "GLFGAIAGFIENGWEGLIDGWYGFRHQNAQ",
                                "GEGTAADYKSTQSAIDQITGKLNRLIEKTN",
                                "QQFELIDNEFTEVEKQIGNVINWTRDSMTE",
                                "VWSYNAELLVAMENQHTIDLADSEMNKLYE",
                                "RVKRQLRENAEEDGTGCFEIFHKCDDDCMA",
                                "SIRNNTYDHSKYREEAMQNRIQINPVKLSS",
                                "GYKDVILWFSFGASCFILLAIAMGLVFICV",
                                "KNGNMRCTICI")
                }
        };

        for (String[] testData : testingData) {

            // must have atleast id, polyprotein, peptide db
            assert (testData.length >= 3);

            ViralProtein protein = new ViralProtein();
            protein.setDefline(testData[0]);
            protein.setProteinID("");
            protein.setSequence(ProteinSequence.of(testData[1]));
            File peptideDB = Paths.get(VigorUtils.getVirusDatabasePath(), testData[2]).toFile();

            List<ProteinSequence> matches = peptideService.findPeptides(protein, peptideDB);
            LOGGER.debug(() -> String.format("peptides:%s", matches.stream().map(String::valueOf).collect(Collectors.joining("\n> ", "\n> ", ""))));

            String[] expected = Arrays.copyOfRange(testData, 3, testData.length);

            assertThat("peptide matches", matches.size(), equalTo(expected.length));
            for (String expectedPeptide : expected) {
                assertThat(matches, hasItem(ProteinSequence.of(expectedPeptide)));
            }
        }
    }

}
